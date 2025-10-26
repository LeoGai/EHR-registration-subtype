rm(list = ls())
library(doParallel)
library("parallel")

generate_data_S5 <- function(num_individuals, time_points, start_time, end_time, feature) {
  data <- data.frame(ID = integer(), Feature = character(), Time = numeric(), Value = numeric())
  
  for (i in 1:num_individuals) {
    
    times <- sort(runif(time_points, min = start_time, max = end_time))
    
    if (i <= 250) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 24 + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else if (i > 250 & i <= 500) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 17 + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else if (i > 500 & i <= 750) {
      base_sequence <- 4*sin(0.9*(times+3))
      individual_sequence <- 28-0.5*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else {
      base_sequence <- 4*sin(0.9*(times+3))
      individual_sequence <- 21-0.5*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    }
    
    individual_data <- data.frame(
      ID = rep(i, length(times)),
      Feature = rep(feature, length(times)),
      Time = times,
      Value = individual_sequence
    )
    
    data <- rbind(data, individual_data)
  }
  
  return(data)
}

apply_time_shifts_S5 <- function(data, num_individuals, timepos_shift, start_time, seed) {
  data_shifted <- data
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(shift_amounts[1:100],rep(0,150),shift_amounts[101:200],rep(0,150),shift_amounts[201:300],rep(0,150),shift_amounts[301:400],rep(0,150) )
  
  for (i in 1:num_individuals) {
    individual_indices <- data_shifted$ID == i
    data_shifted$Time[individual_indices] <- data_shifted$Time[individual_indices] + shift_amounts[i]
    data_shifted <- data_shifted[data_shifted$Time >= start_time & data_shifted$Time <= 17, ]
  }
  
  return(list(data_shifted = data_shifted, shift_amounts = shift_amounts))
}

Tresult <- matrix(nrow = 0, ncol = 5)
colnames(Tresult) <- c("number", "ACC", "LQ1", "MAE", "RUNTIME")

source("./EHR-registration-subtype-main/Simulation_ablation3/Registration_ablation3.R")


cl<-makeCluster(10)
stopCluster(cl)

cls <- makeCluster(40)
registerDoParallel(cls)
clusterExport(cls, "Tresult")
clusterExport(cls, "reg_simu")
clusterExport(cls, "generate_data_S5")
clusterExport(cls, "apply_time_shifts_S5")
Reg_Simu5 <- foreach(seednum = 1:500, .combine = 'rbind', .packages = c("cluster")) %dopar% reg_simu(seednum,generate_data=generate_data_S5,apply_time_shifts=apply_time_shifts_S5)
stopCluster(cls)


calculateStats <- function(matrixData) {
  selectedColumns <- matrixData[, 2:5]
  resultsMatrix <- matrix(nrow = 4, ncol = 4)
  colnames(resultsMatrix) <- colnames(selectedColumns)
  rownames(resultsMatrix) <- c("Mean", "Median", "2.5%", "97.5%")
  for (i in 1:ncol(selectedColumns)) {
    columnData <- selectedColumns[, i]
    resultsMatrix[1, i] <- mean(columnData, na.rm = TRUE) 
    resultsMatrix[2, i] <- median(columnData, na.rm = TRUE) 
    resultsMatrix[3, i] <- quantile(columnData, probs = 0.025, na.rm = TRUE) 
    resultsMatrix[4, i] <- quantile(columnData, probs = 0.975, na.rm = TRUE) 
  }
  resultsMatrix[,4]<-resultsMatrix[,4]/60
  return(resultsMatrix)
}  


Result5 <- calculateStats(Reg_Simu5)

saveRDS(Reg_Simu5, file = paste0("./EHR-registration-subtype-main/Simulation_ablation3/FReg_Simu", 5, ".rds"))

saveRDS(Result5, file = paste0("./EHR-registration-subtype-main/Simulation_ablation3/FMatrix_Simu", 5, ".rds"))

print(Result5)