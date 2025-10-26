rm(list = ls())
library(doParallel)
library(parallel)


reg_simu <- function(seednum) {
  
  library(tidyverse)    
  library(cluster)    
  library(splines)
  library(glmnet)
  source("./EHR-registration-subtype-main/algorithm.R")
  
  seed <- seednum
  set.seed(seed)
  num_individuals <- 1000
  time_points <-28  
  start_time <- 1
  end_time <- 21
  feature <- "A1"
  alpha<-0.95
  k_range <- 2:8
  timepos_option<-c(4,3,2,1)
  timepos_shift<-c(-1,-2,-3,-4)
  
  data <- data.frame(
    ID = integer(),
    Feature = character(),
    Time = numeric(),
    Value = numeric(),
    Scale = numeric()
  )
  
  for (i in 1:num_individuals) {
    t_raw <- sort(runif(time_points, min = start_time, max = end_time))
    
    scale_factor <- 1
    
    if (i <= 300) {
      base_sequence <- 4*sin(0.7*(t_raw+4))+7
      individual_sequence <- 1.5*t_raw + base_sequence + rnorm(length(t_raw), mean = 0, sd =0.8)
    } else if (i > 300 && i <= 550) {
      base_sequence <- 3 * sin(1 * (t_raw + 4))
      individual_sequence <- 2 * t_raw + base_sequence + rnorm(length(t_raw), mean = 0, sd = 0.8)
      scale_factor <- 1
    } else if (i > 550 && i <= 800) {      
      base_sequence <- 3 * sin(1 * (t_raw + 4))
      individual_sequence <- 2 * t_raw + base_sequence + rnorm(length(t_raw), mean = 0, sd = 0.8)
      scale_factor <- 0.7
      
    } else {                                
      base_sequence <- 3 * sin(1 * (t_raw + 4))
      individual_sequence <- 2 * t_raw + base_sequence + rnorm(length(t_raw), mean = 0, sd = 0.8)
      scale_factor <- runif(1, min = 0.7, max = 1)
    }
    
    t_scaled <- t_raw * scale_factor
    
    individual_data <- data.frame(
      ID = rep(i, length(t_scaled)),
      Feature = rep(feature, length(t_scaled)),
      Time = t_scaled,
      Value = individual_sequence,
      Scale = rep(scale_factor, length(t_scaled))
    )
    
    data <- rbind(data, individual_data)
  }
  
  data$GroupCategory <- with(data,
                             ifelse(ID <= 300, "Group 1",
                                    "Group 2")
  )
  
  data <- data[ !(data$ID > 550 & data$Time >= 17 * data$Scale), ]
  
  # Shift time
  data_shifted <- data
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(rep(0,180),shift_amounts[1:120],rep(0,150),shift_amounts[121:220],rep(0,150),shift_amounts[221:320],rep(0,120),shift_amounts[321:400] )
  for (i in 1:num_individuals) {
    individual_indices <- data_shifted$ID == i
    data_shifted$Time[individual_indices] <- data_shifted$Time[individual_indices] + shift_amounts[i]
    data_shifted <- data_shifted[data_shifted$Time >= start_time  & data_shifted$Time <= 17, ]
  }
  
  
  data_recover<-SubtypeAware_Registration(data=data_shifted,alpha=0.95,k_range=2:8,timepos_option=c(4,3,2,1),tau=0.45,tmin=1,tmax=17,nots=c(8,13))
  
  
  library(dtwclust)   
  library(proxy)
  library(mclust)   
  library(aricode) 
  library(clue) 
  
  
  data_recover <- data_recover %>%
    group_by(ID) %>%
    filter(
      sum(Time > 8 & Time <= 13) >= 1
    ) %>% filter(
      sum(Time > 0 & Time <=8) >= 1
    ) %>% 
    ungroup()
  
  
  data_shifted <- data_shifted %>%
    semi_join(data_recover %>% distinct(ID), by = "ID")
  
  data <- data %>%
    semi_join(data_recover %>% distinct(ID), by = "ID")
  
  shifted_data<- data_shifted
  recover_data<- data_recover
  raw_data<-data
  
  series_list <- shifted_data %>%
    arrange(ID, Time) %>%
    group_by(ID) %>%
    group_map(~ as.numeric(.x$Value))
  
  
  ids <- shifted_data %>% distinct(ID) %>% arrange(ID) %>% pull()
  stopifnot(length(series_list) == length(ids))
  keep_idx <- seq_along(series_list)
  X_list  <- series_list[keep_idx]
  ids_sub <- ids[keep_idx]
  n <- length(X_list)
  cat("Will compute pairwise DTW for n =", n, "series\n")
  
  
  W <- round(0.5 * median(sapply(X_list, length)))
  D <- proxy::dist(X_list, method = "dtw_basic", normalize = TRUE, window.size = W)
  k_choose<-2
  k_range <- 2:4
  best_pam <- pam(D, k = k_choose)
  clusters <- best_pam$clustering
  result_df <- data.frame(ID = ids_sub, cluster = clusters)
  
  AKI_data_recover<-recover_data
  AKI_data_recover_21day<-AKI_data_recover
  bspline_coefficients_list_recover_21day<-map_bspline(AKI_data_recover_21day,timepos_option=c(4,3,2,1),tmin=1,tmax=17,nots=c(8,13))
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
  
  best_pam2 <- pam(data_clustering_recover_21day, k = k_choose)
  clusters2 <- best_pam2$clustering
  result_df2 <- data.frame(ID = ids_sub, cluster = clusters2)
  
  AKI_data_raw_21day<-shifted_data
  
  bspline_coefficients_list_21day<-map_bspline(AKI_data_raw_21day,timepos_option=c(4,3,2,1),tmin=1,tmax=17,nots=c(8,13))
  
  coefficients_df0_21day <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())
  
  for (i in 1:length(bspline_coefficients_list_21day)) {
    original_coefficients <- bspline_coefficients_list_21day[[i]]$Original
    coef1 <- original_coefficients["(Intercept)"]
    coef2 <- original_coefficients["bsplines1"]
    coef3 <- original_coefficients["bsplines2"]
    coef4 <- original_coefficients["bsplines3"]
    coef5 <- original_coefficients["bsplines4"]
    coef6 <- original_coefficients["bsplines5"]
    coefficients_df0_21day <- rbind(coefficients_df0_21day, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
  }
  
  data_clustering_raw_21day<- coefficients_df0_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
  
  best_pam3 <- pam(data_clustering_raw_21day, k = k_choose)
  clusters3 <- best_pam3$clustering
  result_df3 <- data.frame(ID = ids_sub, cluster = clusters3)
  
  AKI_data_raw_21day<-raw_data
  
  bspline_coefficients_list_21day<-map_bspline(AKI_data_raw_21day,timepos_option=c(4,3,2,1),tmin=1,tmax=17,nots=c(8,13))
  
  
  coefficients_df0_21day <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())
  
  for (i in 1:length(bspline_coefficients_list_21day)) {
    original_coefficients <- bspline_coefficients_list_21day[[i]]$Original
    coef1 <- original_coefficients["(Intercept)"]
    coef2 <- original_coefficients["bsplines1"]
    coef3 <- original_coefficients["bsplines2"]
    coef4 <- original_coefficients["bsplines3"]
    coef5 <- original_coefficients["bsplines4"]
    coef6 <- original_coefficients["bsplines5"]
    coefficients_df0_21day <- rbind(coefficients_df0_21day, data.frame(ID = i, Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
  }
  
  data_clustering_raw_21day<- coefficients_df0_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
  
  
  best_pam4 <- pam(data_clustering_raw_21day, k = k_choose)
  clusters4 <- best_pam4$clustering
  result_df4 <- data.frame(ID = ids_sub, cluster = clusters4)
  
  ## 1)  ID ture lable
  truth_df <- shifted_data %>%
    distinct(ID, GroupCategory) %>%
    arrange(ID)
  
  ## 
  stopifnot(all(c("ID","cluster") %in% names(result_df2)))
  
  eval_with_truth <- function(pred_df, truth_df) {
    df <- pred_df %>% inner_join(truth_df, by = "ID")
    pred  <- as.integer(as.factor(df$cluster))          
    truth <- as.integer(as.factor(df$GroupCategory))   
    list(pred = pred, truth = truth, df = df)
  }
  
  E1 <- eval_with_truth(result_df,  truth_df)   
  E2 <- eval_with_truth(result_df2, truth_df)   
  
  bestmap_accuracy <- function(pred, truth) {
    tab <- table(pred, truth)
    assignment <- solve_LSAP(tab, maximum = TRUE)
    acc <- sum(tab[cbind(seq_len(nrow(tab)), assignment)]) / length(truth)
    as.numeric(acc)
  }
  
  compute_metrics <- function(pred, truth) {
    list(
      ARI  = adjustedRandIndex(truth, pred),
      AMI  = AMI(truth, pred),
      Acc_bestmap = bestmap_accuracy(pred, truth)
    )
  }
  
  metrics_dtw   <- compute_metrics(E1$pred, E1$truth)
  metrics_spline<- compute_metrics(E2$pred, E2$truth)
  
  
  res_row <- data.frame(
    seed = seed,
    ARI_DTW = metrics_dtw$ARI,
    AMI_DTW = metrics_dtw$AMI,
    ACC_DTW = metrics_dtw$Acc_bestmap,
    ARI_Our = metrics_spline$ARI,
    AMI_Our = metrics_spline$AMI,
    ACC_Our = metrics_spline$Acc_bestmap,
    stringsAsFactors = FALSE
  )
  
  return(res_row)
}

cl <- makeCluster(10)
stopCluster(cl)

cls <- makeCluster(40)
registerDoParallel(cls)


Reg_k2 <- foreach(seednum = 1:500,
                  .combine = 'rbind',
                  .packages = c(
                    "dplyr","dtw","dtwclust","proxy","cluster","mclust",
                    "aricode","clue","glmnet","splines"
                  )) %dopar% {
                    reg_simu(seednum)
                  }
stopCluster(cls)

saveRDS(Reg_k2, file = paste0("./EHR-registration-subtype-main/Comparison_DTW/Reg_k2",".rds"))


## ===============================
library(tidyverse)  
df_long <- Reg_k2 %>%
  pivot_longer(
    cols = -seed,
    names_to = c("metric", "method"),
    names_sep = "_",
    values_to = "value"
  )

summary_stats_k2 <- df_long %>%
  group_by(metric, method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats_k2)

saveRDS(summary_stats_k2, file = paste0("./EHR-registration-subtype-main/Comparison_DTW/summary_stats_k2", ".rds"))
