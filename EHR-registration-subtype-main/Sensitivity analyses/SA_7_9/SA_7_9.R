rm(list = ls())
library(doParallel)
library("parallel")

generate_data_SA79 <- function(num_individuals, time_points, start_time, end_time, feature) {
  data <- data.frame(ID = integer(), Feature = character(), Time = numeric(), Value = numeric())
  for (i in 1:num_individuals) {
    times <- sort(runif(time_points, min = start_time, max = end_time))
    if (i <= 300) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 20 + base_sequence + rnorm(length(times), mean = 0, sd = sqrt(eps[w-6]) )
    } else if (i > 300 & i <= 700) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 17 + base_sequence + rnorm(length(times), mean = 0, sd = sqrt(eps[w-6]) )
    } else {
      base_sequence <- 3*sin(0.9*(times+4))
      individual_sequence <- 16+0.5*times + base_sequence + rnorm(length(times), mean = 0, sd = sqrt(eps[w-6]) )
    }
    individual_data <- data.frame(
      ID = rep(i, length(times)),
      Feature = rep(feature, length(times)),
      Time = times,
      Value = individual_sequence
    )
    data <- rbind(data, individual_data)
  }
  data
}

apply_time_shifts_SA79 <- function(data, num_individuals, timepos_shift, start_time, seed) {
  data_shifted <- data
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(shift_amounts[1:120],rep(0,180),shift_amounts[121:280],rep(0,240),shift_amounts[281:400],rep(0,180) )
  for (i in 1:num_individuals) {
    individual_indices <- data_shifted$ID == i
    data_shifted$Time[individual_indices] <- data_shifted$Time[individual_indices] + shift_amounts[i]
    data_shifted <- data_shifted[data_shifted$Time >= start_time & data_shifted$Time <= 17, ]
  }
  return(list(data_shifted = data_shifted, shift_amounts = shift_amounts))
}


for (w in 7:9) {
  
  eps<-(0.8^2)*c(1.15,1.3,1.45)
  
  Tresult<-matrix(,nrow = 0,ncol=5)
  colnames(Tresult)<-c("number","ACC","LQ1", "MAE", "RUNTIME")
  
  source("./EHR-registration-subtype-main/Sensitivity analyses/Registration_Simulation_SA.R")
  
  
  cl<-makeCluster(10)
  stopCluster(cl)
  
  cls <- makeCluster(40)
  registerDoParallel(cls)
  clusterExport(cls, "Tresult")
  clusterExport(cls, "reg_simu")
  clusterExport(cls, "eps")
  clusterExport(cls, "w")
  clusterExport(cls, "generate_data_SA79")
  clusterExport(cls, "apply_time_shifts_SA79")
  
  saresult <-  foreach(seednum = 1:500, .combine = 'rbind', .packages = c("cluster")) %dopar%  reg_simu(seednum,generate_data=generate_data_SA79,apply_time_shifts=apply_time_shifts_SA79)
  
  
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
  
  
  Result <- calculateStats(saresult)

assign(paste0("Reg_SA", w), saresult)
saveRDS(saresult, file = paste0("./EHR-registration-subtype-main/Sensitivity analyses/SA_7_9/Reg_SA", w, ".rds"))
assign(paste0("Matrix_SA", w), Result)
saveRDS(Result, file = paste0("./EHR-registration-subtype-main/Sensitivity analyses/SA_7_9/Matrix_SA", w, ".rds"))

stopCluster(cls)

}

