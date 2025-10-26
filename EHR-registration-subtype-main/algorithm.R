library(tidyverse)    
library(cluster)    
library(splines)
library(glmnet)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("dist", "stats")
# Function to calculate B-spline coefficients
calculate_bspline_coefficients <- function(data,nots) {
  time_normalized <- data$Time
  bsplines <- bs(time_normalized, knots = nots, degree = 3)
  fit<-glmnet(bsplines, data$Value, alpha = 0, lambda = 0.0000001)
  coefficients <- as.vector(coef(fit))
  names(coefficients) <- c("(Intercept)", paste("bsplines", 1:5, sep=""))
  return(coefficients)
}

############## Initialize B-spline coefficients for each ID################
map_bspline <- function(data,timepos_option,tmin,tmax,nots) {
  bspline_coefficients_list <- list()
  for (id in unique(data$ID)) {
    individual_data <- data[data$ID == id, ]
    original_coefficients <- calculate_bspline_coefficients(individual_data,nots)
    shifted_coefficients <- list()
    for (shift in timepos_option) {
      shifted_data <- individual_data
      shifted_data$Time <- shifted_data$Time + shift
      shifted_data <- shifted_data[shifted_data$Time >= tmin & shifted_data$Time <= tmax, ]
      coefficients <- calculate_bspline_coefficients(shifted_data,nots)
      shifted_coefficients[[paste("Shift", shift)]] <- coefficients
    }
    bspline_coefficients_list[[as.character(id)]] <- list(Original = original_coefficients, Shifted = shifted_coefficients)
  }
  return(bspline_coefficients_list)
}


##################################Main Function##############################

# Input arguments:
#   data : A data frame containing longitudinal observations with columns:
#          - ID: subject identifier
#          - Feature: variable name
#          - Time: observation time 
#          - Value: observed measurement
#          - (Optional) Scale and GroupCategory if available
#
#   alpha : A proportion parameter controlling the selection of the subset for registration, as described in Step 2 of Section 2.
#
#   k_range : Numeric vector specifying the candidate cluster numbers for silhouette-based evaluation (e.g., 2:8).
#
#   timepos_option : Vector defining candidate temporal shift positions (e.g., c(4,3,2,1)).
#
#   tau : Numeric, predefined threshold for the average silhouette coefficient (e.g., 0.45). 
#         As described in Step 2 of Section 2.
#
#   tmin, tmax : Lower and upper bounds of the observed time range.
#
#   nots : Numeric vector specifying knot positions for temporal spline fitting (e.g., c(8,13)).


# Output:
#   Returns a data frame of the same structure as the input (data), with recovered (registered) temporal positions.




SubtypeAware_Registration<-function(data,alpha,k_range,timepos_option,tau,tmin,tmax,nots) {
  #,seeds
  #set.seed(seeds)
  bspline_coefficients_list<-map_bspline(data,timepos_option,tmin,tmax,nots)
  
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
    km_res <- pam(data_clustering, k = k)
    silhouette_res <- silhouette(km_res$clustering, dist(data_clustering))
    silhouette_widths[k - min(k_range) + 1] <- mean(silhouette_res[, "sil_width"])
  }
  
  optimal_k <- k_range[which.max(silhouette_widths)]
  print(silhouette_widths)
  final_clustering_result <- pam(data_clustering, k = optimal_k)
  coefficients_df0$Cluster <- final_clustering_result$clustering
  
  selected_points_list <- list()
  
  for (k in 1:optimal_k) {
    cluster_data <- filter(coefficients_df0, Cluster == k)
    distances <- sqrt(rowSums((cluster_data[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - final_clustering_result$medoids[k, ])^2))
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
      local_center <- colMeans(selected_points[,c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
      
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
  
  #################################Condition of iteration###########################################
  
  if(max(silhouette_widths)>tau){
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
      km_res <- pam(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")],  k = k)
      silhouette_res <- silhouette(km_res$clustering, dist(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]))
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
    if (iteration_count > 50) {
      break
    }
    previous_k <- optimal_k
    s_old<-silhouette_widths
    final_clustering_result <- pam(data_clustering_adj[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")],k = optimal_k)
    data_clustering_adj$Cluster <- final_clustering_result$clustering
    selected_points_list <- list()
    
    for (k in 1:optimal_k) {
      cluster_data <- filter(data_clustering_adj, Cluster == k)
      distances <- sqrt(rowSums((cluster_data[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")] - final_clustering_result$medoids[k, ])^2))
      threshold_distance <- quantile(distances,alpha)
      selected_points <- cluster_data[distances <= threshold_distance, ]
      selected_points_list[[k]] <- selected_points
    }
    
    selected_points_centers <- lapply(selected_points_list, function(cluster_points) {
      colMeans(cluster_points[,c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
    })
    
    optimal_states_results <- data.frame(ID = integer(), Cluster = integer(), OptimalState = character(), MinDistance = numeric(), stringsAsFactors = FALSE)
    
    
    for (k in 1:optimal_k) {
      repeat {
        
        any_update <- FALSE
        selected_points <- selected_points_list[[k]]
        local_center <- colMeans(selected_points[,c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")], na.rm = TRUE)
        
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
  
  ############################## Step to classify remaining points########################################3
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
      
      center_coords <- cluster_centers[center_idx,  c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5","Coef6")]
      distances <- apply(remaining_infor[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5","Coef6")], 1, function(coefs) {
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
  AKI_data_recover  <- data
  shift_recover <- data_optimal_states_result_final$ShiftNumber
  
  for (i in 1:length(bspline_coefficients_list)) {
    
    individual_indices <-AKI_data_recover$ID == unique(AKI_data_recover$ID)[i]
    AKI_data_recover$Time[individual_indices] <- AKI_data_recover$Time[individual_indices] + shift_recover[i]
    AKI_data_recover <-AKI_data_recover[AKI_data_recover$Time >= tmin & AKI_data_recover$Time <= tmax, ]
  }
  
  return(AKI_data_recover)
  
  
}

####Example##########
#AKI_data_recover<-SubtypeAware_Registration(data=AKI_data_raw_filtered,alpha=0.95,k_range=2:8,timepos_option=c(4,3,2,1),tau=0.45,tmin=0,tmax=21,nots=c(7,14),seeds=1111)



