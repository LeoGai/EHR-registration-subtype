library(tidyverse)
library(cluster)
library(gridExtra)
library(splines)
library(glmnet)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
set.seed(6666)
alpha<-0.95
k_range <- 2:8
timepos_option<-c(4,3,2,1)

AKI_diag<-read.csv("./EHR-registration-subtype-main/Real_data/diagnose_AKI_May_iv.csv")
uni_id<-unique(AKI_diag$ICU_ID)
AKI_lab<-read.csv("./EHR-registration-subtype-main/Real_data/lab_AKI_May_iv.csv")
AKI_base<-read.csv("./EHR-registration-subtype-main/Real_data/AKI_raw_base_0430.csv")
endstage_index <- grepl("End stage renal disease", AKI_diag$DIAGNOSE, ignore.case = TRUE)
# drop End stage renal disease
endstage_output <- AKI_diag$ICU_ID[endstage_index ]
drop_endstage_id<-unique(endstage_output)
uni_id_01 <- uni_id[!uni_id %in% drop_endstage_id]
# drop burn
burn_index <- grepl("(?i)\\b(?!heartburn)\\w*burn\\b", AKI_diag$DIAGNOSE, perl = TRUE)
burn_output <- AKI_diag$ICU_ID[burn_index ]
uni_id_02 <-uni_id_01 [!uni_id_01  %in% unique(burn_output)]
# drop renal dialysis
rd_index <- grepl("renal dialysis", AKI_diag$DIAGNOSE, ignore.case = TRUE)
rd_output <- AKI_diag$ICU_ID[rd_index ]
uni_id_03 <-uni_id_02[!uni_id_02  %in% unique(rd_output)]
#### drop estimated glomerular filtration rate <15
AKI_eGFR<-read.csv("./EHR-registration-subtype-main/Real_data/lab_eGFR_AKI_Apr_iv_new.csv")
icu_ids_with_first_value_below_15 <- AKI_eGFR %>%
  group_by(ICU_ID) %>%
  filter(row_number() == 1 & VALUE < 15) %>%
  ungroup() %>%
  select(ICU_ID)

low_gfr_id <- icu_ids_with_first_value_below_15$ICU_ID


#########################################################
drop_icuid<- setdiff(uni_id, uni_id_03)

AKI_lab_SC  <- AKI_lab %>%
  filter(FEATURE_NAME == "Creatinine", !(ICU_ID %in% drop_icuid))

AKI_lab_SC  <- AKI_lab_SC  %>%
  filter(FEATURE_NAME == "Creatinine", !(ICU_ID %in% low_gfr_id))


AKI_data_raw<-AKI_lab_SC[,c(5,2,3,4,1)]
colnames(AKI_data_raw)<-c("ID","Feature","Time","Value","SUBJECT_ID")
length(unique(AKI_data_raw[,1]))
AKI_data_raw[,3]<-AKI_data_raw[,3]/(60*24)
AKI_data_raw_filtered <- AKI_data_raw %>%
  filter(Time <= 21)

AKI_data_raw_filtered <- AKI_data_raw_filtered %>%
  group_by(ID) %>%
  mutate(count = n()) %>%
  filter(count >= 15) %>%
  ungroup()

AKI_data_raw_filtered$BASEVALUE <- AKI_base$Base_SCr

AKI_data_raw_filtered <- AKI_data_raw_filtered %>%
  group_by(ID) %>%
  filter(
    sum(Time >= 0 & Time <= 7) >= 3 &   
      sum(Time > 7 & Time <= 14) >= 3 &   
      sum(Time > 14 & Time <= 21) >= 3    
  ) %>%
  ungroup()

AKI_data_raw_filtered$count <- NULL 
AKI_data_raw_filtered<-AKI_data_raw_filtered[!is.na(AKI_data_raw_filtered[, 4]), ]


time1<-Sys.time()

# Function to calculate B-spline coefficients
calculate_bspline_coefficients <- function(data) {
  time_normalized <- data$Time
  bsplines <- bs(time_normalized, knots = c(7,14), degree = 3)
  fit<-glmnet(bsplines, data$Value, alpha = 0, lambda = 0.0000001)
  coefficients <- as.vector(coef(fit))
  names(coefficients) <- c("(Intercept)", paste("bsplines", 1:5, sep=""))
  return(coefficients)
}

# Initialize B-spline coefficients for each ID
map_bspline <- function(data) {
  bspline_coefficients_list <- list()
  for (id in unique(data$ID)) {
    individual_data <- data[data$ID == id, ]
    original_coefficients <- calculate_bspline_coefficients(individual_data)
    shifted_coefficients <- list()
    for (shift in timepos_option) {
      shifted_data <- individual_data
      shifted_data$Time <- shifted_data$Time + shift
      shifted_data <- shifted_data[shifted_data$Time >= 0 & shifted_data$Time <= 21, ]
      coefficients <- calculate_bspline_coefficients(shifted_data)
      shifted_coefficients[[paste("Shift", shift)]] <- coefficients
    }
    bspline_coefficients_list[[as.character(id)]] <- list(Original = original_coefficients, Shifted = shifted_coefficients)
  }
  return(bspline_coefficients_list)
}

bspline_coefficients_list<-map_bspline(AKI_data_raw_filtered)

coefficients_df <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())

for (i in 1:length(bspline_coefficients_list)) {
  original_coefficients <- bspline_coefficients_list[[i]]$Original
  coef1 <- original_coefficients["(Intercept)"]
  coef2 <- original_coefficients["bsplines1"]
  coef3 <- original_coefficients["bsplines2"]
  coef4 <- original_coefficients["bsplines3"]
  coef5 <- original_coefficients["bsplines4"]
  coef6 <- original_coefficients["bsplines5"]
  coefficients_df <- rbind(coefficients_df, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
}

coefficients_df0 <- coefficients_df


# Shifted data
coefficients_dfp1 <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())
coefficients_dfp2 <- coefficients_dfp1
coefficients_dfp3 <- coefficients_dfp1
coefficients_dfp4 <- coefficients_dfp1

for (i in 1:length(bspline_coefficients_list)) {
  shifted_coefficients_1 <- bspline_coefficients_list[[i]]$Shifted$`Shift 1`
  shifted_coefficients_2 <- bspline_coefficients_list[[i]]$Shifted$`Shift 2`
  shifted_coefficients_3 <- bspline_coefficients_list[[i]]$Shifted$`Shift 3`
  shifted_coefficients_4 <- bspline_coefficients_list[[i]]$Shifted$`Shift 4`
  
  for (shifted_coefficients in list(shifted_coefficients_1, shifted_coefficients_2, shifted_coefficients_3, shifted_coefficients_4)) {
    coef1 <- shifted_coefficients["(Intercept)"]
    coef2 <- shifted_coefficients["bsplines1"]
    coef3 <- shifted_coefficients["bsplines2"]
    coef4 <- shifted_coefficients["bsplines3"]
    coef5 <- shifted_coefficients["bsplines4"]
    coef6 <- shifted_coefficients["bsplines5"]
    coeff_dfp <- data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6)
    
    if (identical(shifted_coefficients, shifted_coefficients_1)) {
      coefficients_dfp1 <- rbind(coefficients_dfp1, coeff_dfp)
    } else if (identical(shifted_coefficients, shifted_coefficients_2)) {
      coefficients_dfp2 <- rbind(coefficients_dfp2, coeff_dfp)
    } else if (identical(shifted_coefficients, shifted_coefficients_3)) {
      coefficients_dfp3 <- rbind(coefficients_dfp3, coeff_dfp)
    } else if (identical(shifted_coefficients, shifted_coefficients_4)) {
      coefficients_dfp4 <- rbind(coefficients_dfp4, coeff_dfp)
    }
  }
}


coefficients_df0$State <- "Original"
coefficients_dfp1$State <- "Shift 1"
coefficients_dfp2$State <- "Shift 2"
coefficients_dfp3$State <- "Shift 3"
coefficients_dfp4$State <- "Shift 4"
all_coefficients_df <- rbind(coefficients_df0, coefficients_dfp1, coefficients_dfp2, coefficients_dfp3, coefficients_dfp4)


data_clustering <- coefficients_df0[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
silhouette_widths <- numeric(length(k_range))

for (k in k_range) {
  km_res <- kmeans(data_clustering, centers = k, nstart = 25)
  silhouette_res <- silhouette(km_res$cluster, dist(data_clustering))
  silhouette_widths[k - min(k_range) + 1] <- mean(silhouette_res[, "sil_width"])
}

optimal_k <- k_range[which.max(silhouette_widths)]
print(silhouette_widths)
final_clustering_result <- kmeans(data_clustering, centers = optimal_k, nstart = 50)
coefficients_df0$Cluster <- final_clustering_result$cluster

selected_points_list <- list()

for (k in 1:optimal_k) {
  cluster_data <- filter(coefficients_df0, Cluster == k)
  distances <- sqrt(rowSums((cluster_data[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - final_clustering_result$centers[k, ])^2))
  threshold_distance <- quantile(distances, alpha)
  selected_points <- cluster_data[distances <= threshold_distance, ]
  selected_points_list[[k]] <- selected_points
}


selected_points_centers <- lapply(selected_points_list, function(cluster_points) {
  colMeans(cluster_points[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
})


optimal_states_results <- data.frame(ID = integer(), Cluster = integer(), OptimalState = character(), MinDistance = numeric(), stringsAsFactors = FALSE)

for (k in 1:optimal_k) {
  repeat {
    any_update <- FALSE
    selected_points <- selected_points_list[[k]]
    local_center <- colMeans(selected_points[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
    
    for (i in 1:nrow(selected_points)) {
      point_id <- selected_points$ID[i]
      point_data_all_states <- all_coefficients_df[all_coefficients_df$ID == point_id, ]
      distances <- sapply(1:nrow(point_data_all_states), function(row) {
        sum((point_data_all_states[row, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - local_center)^2)
      })
      
      min_distance_index <- which.min(distances)
      optimal_state <- point_data_all_states$State[min_distance_index]
      min_distance <- distances[min_distance_index]
      current_state <- selected_points$State[i]
      if (optimal_state != current_state) {
        
        selected_points$State[i] <- optimal_state
        selected_points_list[[k]] <- selected_points  
        any_update <- TRUE  
        optimal_states_results <- rbind(optimal_states_results, data.frame(ID = point_id, Cluster = k, OptimalState = optimal_state, MinDistance = min_distance))
      }
      if (optimal_state == current_state){
        optimal_states_results <- rbind(optimal_states_results, data.frame(ID = point_id, Cluster = k, OptimalState = current_state, MinDistance = min_distance))
      }
    }
    
    if (!any_update) {
      break
    }
  }
}



optimal_states_results <- distinct(optimal_states_results, ID, .keep_all = TRUE)

optimal_states_results <- optimal_states_results %>%
  mutate(ShiftNumber = case_when(
    OptimalState == "Original" ~ 0,
    OptimalState == "Shift -3" ~ -3,
    OptimalState == "Shift -2" ~ -2,
    OptimalState == "Shift -1" ~ -1,
    OptimalState == "Shift 1" ~ 1,
    OptimalState == "Shift 2" ~ 2,
    OptimalState == "Shift 3" ~ 3,
    OptimalState == "Shift 4" ~ 4,
    OptimalState == "Shift 5" ~ 5,
    OptimalState == "Shift 6" ~ 6,
    OptimalState == "Shift 7" ~ 7,
    TRUE ~ NA_integer_
  ))



all_coefficients_df$ID <- as.integer(all_coefficients_df$ID)
optimal_states_results$ID <- as.integer(optimal_states_results$ID)
all_coefficients_df_adj <- left_join(all_coefficients_df, optimal_states_results[, c("ID", "OptimalState")], by = "ID")
all_coefficients_df_adj <- all_coefficients_df_adj %>%
  filter(State == OptimalState)

data_clustering_adj<-coefficients_df0
data_clustering_adj[all_coefficients_df_adj$ID,1:8]<-all_coefficients_df_adj[,1:8]
#data_clustering_adj<-data_clustering_adj[,-9]


data_optimal_states_result_one<-data_clustering_adj %>%
  mutate(ShiftNumber = case_when(
    State == "Original" ~ 0,
    State == "Shift -3" ~ -3,
    State == "Shift -2" ~ -2,
    State == "Shift -1" ~ -1,
    State == "Shift 1" ~ 1,
    State == "Shift 2" ~ 2,
    State == "Shift 3" ~ 3,
    State == "Shift 4" ~ 4,
    State == "Shift 5" ~ 5,
    State == "Shift 6" ~ 6,
    State == "Shift 7" ~ 7,
    TRUE ~ NA_integer_
  ))

#################################Condition of iteration##################

# Iterative optimization
if(max(silhouette_widths)>0.55){
  max_iterations<-1
}else{
  max_iterations<-100
}

iteration_count <- 0
previous_k <- optimal_k
s_old<-silhouette_widths

repeat {
  iteration_count <- iteration_count + 1
  silhouette_widths <- numeric(length(k_range))
  
  for (k in k_range) {
    set.seed(6)
    km_res <- kmeans(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], centers = k, nstart = 25)
    silhouette_res <- silhouette(km_res$cluster, dist(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]))
    silhouette_widths[k - min(k_range) + 1] <- mean(silhouette_res[, "sil_width"])
  }
  
  
  print(silhouette_widths)
  optimal_k <- k_range[which.max(silhouette_widths)]
  
  
  
  if (optimal_k == previous_k && all(sort(silhouette_widths, decreasing = TRUE)[1:2] <= sort(s_old, decreasing = TRUE)[1:2])) {
    break
  }
  if (max_iterations == 1) {
    break
  }
  if (iteration_count > 20) {
    break
  }
  previous_k <- optimal_k
  s_old<-silhouette_widths
  final_clustering_result <- kmeans(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], centers = optimal_k, nstart = 50)
  data_clustering_adj$Cluster <- final_clustering_result$cluster
  selected_points_list <- list()
  
  for (k in 1:optimal_k) {
    cluster_data <- filter(data_clustering_adj, Cluster == k)
    distances <- sqrt(rowSums((cluster_data[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - final_clustering_result$centers[k, ])^2))
    threshold_distance <- quantile(distances,alpha)
    selected_points <- cluster_data[distances <= threshold_distance, ]
    selected_points_list[[k]] <- selected_points
  }
  
  selected_points_centers <- lapply(selected_points_list, function(cluster_points) {
    colMeans(cluster_points[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
  })
  
  optimal_states_results <- data.frame(ID = integer(), Cluster = integer(), OptimalState = character(), MinDistance = numeric(), stringsAsFactors = FALSE)
  
  for (k in 1:optimal_k) {
    repeat {
      
      any_update <- FALSE
      selected_points <- selected_points_list[[k]]
      local_center <- colMeans(selected_points[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
      
      for (i in 1:nrow(selected_points)) {
        point_id <- selected_points$ID[i]
        point_data_all_states <- all_coefficients_df[all_coefficients_df$ID == point_id, ]
        distances <- sapply(1:nrow(point_data_all_states), function(row) {
          sum((point_data_all_states[row, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - local_center)^2)
        })
        
        min_distance_index <- which.min(distances)
        optimal_state <- point_data_all_states$State[min_distance_index]
        min_distance <- distances[min_distance_index]
        current_state <- selected_points$State[i]
        if (optimal_state != current_state) {
          selected_points$State[i] <- optimal_state
          selected_points_list[[k]] <- selected_points  
          any_update <- TRUE  
          optimal_states_results <- rbind(optimal_states_results, data.frame(ID = point_id, Cluster = k, OptimalState = optimal_state, MinDistance = min_distance))
        }
        if (optimal_state == current_state){
          optimal_states_results <- rbind(optimal_states_results, data.frame(ID = point_id, Cluster = k, OptimalState = current_state, MinDistance = min_distance))
        }
      }
      
      if (!any_update) {
        break
      }
    }
    
    
    
  }
  optimal_states_results <- distinct(optimal_states_results, ID, .keep_all = TRUE)
  all_coefficients_df$ID <- as.integer(all_coefficients_df$ID)
  optimal_states_results$ID <- as.integer(optimal_states_results$ID)
  all_coefficients_df_adj <- left_join(all_coefficients_df, optimal_states_results[, c("ID", "OptimalState")], by = "ID")
  all_coefficients_df_adj <- all_coefficients_df_adj %>%
    filter(State == OptimalState)
  data_clustering_adj[all_coefficients_df_adj$ID,1:8]<-all_coefficients_df_adj[,1:8]
  
}


optimal_states_results <- optimal_states_results %>%
  mutate(ShiftNumber = case_when(
    OptimalState == "Original" ~ 0,
    OptimalState == "Shift -3" ~ -3,
    OptimalState == "Shift -2" ~ -2,
    OptimalState == "Shift -1" ~ -1,
    OptimalState == "Shift 1" ~ 1,
    OptimalState == "Shift 2" ~ 2,
    OptimalState == "Shift 3" ~ 3,
    OptimalState == "Shift 4" ~ 4,
    OptimalState == "Shift 5" ~ 5,
    OptimalState == "Shift 6" ~ 6,
    OptimalState == "Shift 7" ~ 7,
    TRUE ~ NA_integer_
  ))

data_optimal_states_result<-data_clustering_adj %>%
  mutate(ShiftNumber = case_when(
    State == "Original" ~ 0,
    State == "Shift -3" ~ -3,
    State == "Shift -2" ~ -2,
    State == "Shift -1" ~ -1,
    State == "Shift 1" ~ 1,
    State == "Shift 2" ~ 2,
    State == "Shift 3" ~ 3,
    State == "Shift 4" ~ 4,
    State == "Shift 5" ~ 5,
    State == "Shift 6" ~ 6,
    State == "Shift 7" ~ 7,
    TRUE ~ NA_integer_
  ))

#############################Step to classify remaining points########################################3
remaining_points <- data_optimal_states_result %>% 
  filter(!(ID %in% optimal_states_results$ID))

cluster_centers <- data_optimal_states_result %>% 
  filter(ID %in% optimal_states_results$ID) %>%
  group_by(Cluster) %>%
  summarise(Coef1 = mean(Coef1), Coef2 = mean(Coef2), Coef3 = mean(Coef3), Coef4 = mean(Coef4), Coef5 = mean(Coef5), Coef6 = mean(Coef6))

Remain_info <- data.frame(ID = integer(), Cluster = integer(), MinDistance = numeric(), BestState = character())

for (j in remaining_points$ID) {
  remaining_infor<-all_coefficients_df[all_coefficients_df$ID == j, ]
  best_state_info <- data.frame(ID = integer(), Cluster = integer(), MinDistance = numeric(), BestState = character())
  
  
  for (center_idx in 1:nrow(cluster_centers)) {
    
    center_coords <- cluster_centers[center_idx, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
    distances <- apply(remaining_infor[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], 1, function(coefs) {
      sum((coefs -  center_coords)^2)
    })
    distances<- na.omit(distances)
    min_distance <- min(distances)
    min_distance_idx <- which.min(distances)
    best_state <- remaining_infor$State[min_distance_idx]
    
    if (nrow(best_state_info) == 0 || min_distance < best_state_info$MinDistance) {
      best_state_info <- data.frame(ID = j, Cluster = center_idx, MinDistance = min_distance, BestState = best_state)
    }
    
  }
  Remain_info<-rbind(Remain_info,best_state_info)
}

data_optimal_states_result_final <- data_optimal_states_result
data_optimal_states_result_final[Remain_info$ID,"State" ]<-Remain_info$BestState


data_optimal_states_result_final<-data_optimal_states_result_final %>%
  mutate(ShiftNumber = case_when(
    State == "Original" ~ 0,
    State == "Shift -3" ~ -3,
    State == "Shift -2" ~ -2,
    State == "Shift -1" ~ -1,
    State == "Shift 1" ~ 1,
    State == "Shift 2" ~ 2,
    State == "Shift 3" ~ 3,
    State == "Shift 4" ~ 4,
    State == "Shift 5" ~ 5,
    State == "Shift 6" ~ 6,
    State == "Shift 7" ~ 7,
    TRUE ~ NA_integer_
  ))

###########
AKI_data_recover2  <- AKI_data_raw_filtered
shift_recover <- data_optimal_states_result_final$ShiftNumber

for (i in 1:length(bspline_coefficients_list)) {
  
  individual_indices <-AKI_data_recover2$ID == unique(AKI_data_recover2$ID)[i]
  AKI_data_recover2$Time[individual_indices] <- AKI_data_recover2$Time[individual_indices] + shift_recover[i]
  AKI_data_recover2 <-AKI_data_recover2[AKI_data_recover2$Time >= 0 & AKI_data_recover2$Time <= 21, ]
}


run_time_rda <- as.numeric(difftime(Sys.time(), time1, units = "secs"))
print(run_time_rda)

saveRDS(AKI_data_recover2,"./EHR-registration-subtype-main/Real_data/AKI_DA_KM_RE/AKI_data_recover2.rds")


write.csv(AKI_data_recover2, file = "./EHR-registration-subtype-main/Real_data/AKI_DA_KM_RE/AKI_data_recover2.csv", row.names = FALSE)


###################Unsupervised learning task#######################

#######Under 21 days##################
AKI_data_raw_21day<-AKI_data_raw_filtered
bspline_coefficients_list_21day<-map_bspline(AKI_data_raw_21day)


coefficients_df0_21day <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())

for (i in 1:length(bspline_coefficients_list_21day)) {
  original_coefficients <- bspline_coefficients_list_21day[[i]]$Original
  coef1 <- original_coefficients["(Intercept)"]
  coef2 <- original_coefficients["bsplines1"]
  coef3 <- original_coefficients["bsplines2"]
  coef4 <- original_coefficients["bsplines3"]
  coef5 <- original_coefficients["bsplines4"]
  coef6 <- original_coefficients["bsplines5"]
  coefficients_df0_21day <- rbind(coefficients_df0_21day, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))#, Coef6 = coef6))
}

data_clustering_raw_21day<- coefficients_df0_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]

for (j in 2:8) {
  set.seed(6666)
  km_res_raw_21day <-kmeans(data_clustering_raw_21day, centers = j, nstart = 25)
  silhouette_res_raw_21day <- silhouette(km_res_raw_21day$cluster, dist(data_clustering_raw_21day))
  print("================raw====================")
  print(round(mean(silhouette_res_raw_21day[, "sil_width"]),3))
  
}

AKI_data_recover_21day<-AKI_data_recover2
bspline_coefficients_list_recover_21day<-map_bspline(AKI_data_recover_21day)
coefficients_df0_recover_21day <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())

for (i in 1:length(bspline_coefficients_list_recover_21day)) {
  original_coefficients <- bspline_coefficients_list_recover_21day[[i]]$Original
  coef1 <- original_coefficients["(Intercept)"]
  coef2 <- original_coefficients["bsplines1"]
  coef3 <- original_coefficients["bsplines2"]
  coef4 <- original_coefficients["bsplines3"]
  coef5 <- original_coefficients["bsplines4"]
  coef6 <- original_coefficients["bsplines5"]
  coefficients_df0_recover_21day <- rbind(coefficients_df0_recover_21day, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
}

data_clustering_recover_21day <- coefficients_df0_recover_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]

for (j in 2:8) {
  set.seed(6666)
  km_res_recover_21day <- kmeans(data_clustering_recover_21day, centers = j, nstart = 25)
  silhouette_res_recover_21day <- silhouette(km_res_recover_21day$cluster, dist(data_clustering_recover_21day))
  print("==============recover======================")
  print(round(mean(silhouette_res_recover_21day[, "sil_width"]),3))
  
}



#######Under 7 days##################

calculate_bspline_coefficients_7d <- function(data) {
  time_normalized <-data$Time
  bsplines <- bs(time_normalized, knots = c(2,4.5),degree = 3)
  fit<-glmnet(bsplines, data$Value, alpha = 0, lambda = 0.0000001)
  coefficients <- as.vector(coef(fit))
  names(coefficients) <- c("(Intercept)", paste("bsplines", 1:5, sep=""))
  return(coefficients)
}

timepos_option<-c(1)
map_bspline_7d <- function(data) {
  
  bspline_coefficients_list <- list()
  for (id in unique(data$ID)) {
    individual_data <- data[data$ID == id, ]
    original_coefficients <- calculate_bspline_coefficients_7d(individual_data)
    bspline_coefficients_list[[as.character(id)]] <- list(Original = original_coefficients)
  }
  
  return(bspline_coefficients_list)
  
}


#############################################Clinical setting task NEW#############################################

tar_day<-10



AKI_recover_SC<-AKI_data_recover2[,c(1:4,6)]
colnames(AKI_recover_SC)<-c("ICU_ID", "FEATURE_NAME", "RECORD_MIN", "VALUE","BASEVALUE")

AKI_recover_SC <- na.omit(AKI_recover_SC)

AKI_lab_SC_adj_recover <- AKI_recover_SC %>% 
  arrange(ICU_ID, RECORD_MIN) %>%  
  group_by(ICU_ID) %>%  
  filter(RECORD_MIN <= 7) %>%
  mutate(BASELINE_VALUE = BASEVALUE,  
         BASELINE_TIME = 0, 
         VALUE_RATIO = VALUE / BASELINE_VALUE) %>%  
  mutate(
    MAX_VALUE_2d = ifelse(RECORD_MIN <= 2, max(VALUE[RECORD_MIN <= 2], na.rm = TRUE), NA),
    Increase_VALUE_2d=ifelse(RECORD_MIN <= 2, (max(VALUE[RECORD_MIN <= 2], na.rm = TRUE)-BASELINE_VALUE), NA),
    MAX_VALUE_RATIO_7d = ifelse(RECORD_MIN <= 7, max(VALUE_RATIO, na.rm = TRUE), NA)) %>%  
  mutate(SC_CATEGORY = case_when( 
    BASELINE_TIME <= 2 & (MAX_VALUE_2d >= 4.0) ~ "S4",
    BASELINE_TIME <= 2 & (MAX_VALUE_RATIO_7d >= 3) ~ "S3",
    BASELINE_TIME <= 2 & MAX_VALUE_RATIO_7d >= 2 & MAX_VALUE_RATIO_7d < 3 ~ "S2", 
    BASELINE_TIME <= 2 & MAX_VALUE_RATIO_7d >= 1.5 & MAX_VALUE_RATIO_7d < 2 ~ "S1",
    BASELINE_TIME <= 2 & (Increase_VALUE_2d >= 0.3) ~ "S0",
    TRUE ~ "unknown" 
  )) %>%
  summarise(SC_CATEGORY = first(SC_CATEGORY)) 


names(AKI_lab_SC_adj_recover) <- c("ID", "SC_CATEGORY")


AKI_data_recover_Clinic<-AKI_data_recover2 %>%
  inner_join(AKI_lab_SC_adj_recover, by = "ID")

AKI_data_recover_Clinic <- AKI_data_recover_Clinic  %>%
  filter(SC_CATEGORY != "unknown")


AKI_data_recover_Clinic <- AKI_data_recover_Clinic  %>%
  filter(SC_CATEGORY != "S0")


AKI_data_recover_Clinic_7day<-AKI_data_recover_Clinic%>%
  filter(Time <= tar_day)


bspline_coefficients_list_Clinic_7day_recover<-map_bspline_7d(AKI_data_recover_Clinic_7day)


coefficients_df0_Clinic_7day_recover <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())

for (i in 1:length(bspline_coefficients_list_Clinic_7day_recover)) {
  original_coefficients <- bspline_coefficients_list_Clinic_7day_recover[[i]]$Original
  coef1 <- original_coefficients["(Intercept)"]
  coef2 <- original_coefficients["bsplines1"]
  coef3 <- original_coefficients["bsplines2"]
  coef4 <- original_coefficients["bsplines3"]
  coef5 <- original_coefficients["bsplines4"]
  coef6 <- original_coefficients["bsplines5"]
  coefficients_df0_Clinic_7day_recover <- rbind(coefficients_df0_Clinic_7day_recover, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
}

data_clustering_Clinic_7day_recover <- coefficients_df0_Clinic_7day_recover[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]



data_clustering_Clinic_7day_recover[,7]<-unique(AKI_data_recover_Clinic_7day$ID)  
colnames(data_clustering_Clinic_7day_recover)[7]<-"ID"




comb_group<-matrix(c(1, 2, 3, 4, 1, 2, 3, 3,1,1,2,3,1,1,2,2,1,1,1,2), nrow = 5, byrow = TRUE)


for (g in 1:5) {
  
  data_clustering_Clinic_7day_sc_recover<- data_clustering_Clinic_7day_recover %>%
    inner_join(AKI_lab_SC_adj_recover, by = "ID")  
  
  data_clustering_Clinic_7day_sc_recover$SC_CATEGORY <- recode(data_clustering_Clinic_7day_sc_recover$SC_CATEGORY,
                                                               
                                                               "S1" = comb_group[g,1],
                                                               "S2" = comb_group[g,2],
                                                               "S3" = comb_group[g,3],
                                                               "S4" = comb_group[g,4])
  
  silhouette_res_Clinic_recover <- silhouette(data_clustering_Clinic_7day_sc_recover$SC_CATEGORY , dist(data_clustering_Clinic_7day_sc_recover[,1:6]))
  cat("Recovered Data Group",g,":")
  cat(as.numeric(round(mean(silhouette_res_Clinic_recover[, "sil_width"]),3)), "\n")
  
  
}



