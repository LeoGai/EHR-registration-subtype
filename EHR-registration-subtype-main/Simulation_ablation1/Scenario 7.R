rm(list = ls())
library(doParallel)
library("parallel")

generate_data_S7 <- function(num_individuals, time_points, start_time, end_time, feature) {
  data <- data.frame(ID = integer(), Feature = character(), Time = numeric(), Value = numeric())
  
  for (i in 1:num_individuals) {
    times <- sort(runif(time_points, min = start_time, max = end_time))
    if (i <= 500) {
      base_sequence <- 2*sin(0.6*(times+4))
      individual_sequence <- 17 + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else if (i > 500 & i <= 750) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 2*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    }else {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 2*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
      times<-times*0.7
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

apply_time_shifts_S7 <- function(data, num_individuals, timepos_shift, start_time, seed) {
  data_shifted <- data
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(rep(0,300),shift_amounts[1:200],rep(0,150),shift_amounts[201:300],rep(0,150),shift_amounts[301:400] )
  for (i in 1:num_individuals) {
    individual_indices <- data_shifted$ID == i
    data_shifted$Time[individual_indices] <- data_shifted$Time[individual_indices] + shift_amounts[i]
    data_shifted <- data_shifted[data_shifted$Time >= start_time & data_shifted$Time <= 17, ]
  }
  
  return(list(data_shifted = data_shifted, shift_amounts = shift_amounts))
}

Tresult<-matrix(,nrow = 0,ncol=5)
colnames(Tresult)<-c("number","ACC","LQ1", "MAE", "RUNTIME")

source("./EHR-registration-subtype-main/Simulation_ablation1/Registration_ablation1.R")

cl<-makeCluster(10)
stopCluster(cl)

cls <- makeCluster(40)
registerDoParallel(cls)
clusterExport(cls, "Tresult")
clusterExport(cls, "reg_simu")

Reg_Simu7 <- foreach(seednum = 1:500, .combine = 'rbind', .packages = c("cluster")) %dopar% reg_simu(seednum,generate_data=generate_data_S7,apply_time_shifts=apply_time_shifts_S7)
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


Result7 <- calculateStats(Reg_Simu7)

saveRDS(Reg_Simu7, file = paste0("./EHR-registration-subtype-main/Simulation_ablation1/FReg_Simu", 7, ".rds"))

saveRDS(Result7, file = paste0("./EHR-registration-subtype-main/Simulation_ablation1/FMatrix_Simu", 7, ".rds"))


print(Result7)