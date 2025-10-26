rm(list = ls())
library(tidyverse)    
library(cluster)    
library(splines)
library(glmnet)
library(ggplot2)
library(gridExtra)
library(patchwork)
  
  # Function to calculate B-spline coefficients
  calculate_bspline_coefficients <- function(data) {
    time_normalized <- data$Time
    bsplines <- bs(time_normalized, knots = c(8,13), degree = 3)
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
        shifted_data <- shifted_data[shifted_data$Time >= start_time & shifted_data$Time <=17, ]
        coefficients <- calculate_bspline_coefficients(shifted_data)
        shifted_coefficients[[paste("Shift", shift)]] <- coefficients
      }
      bspline_coefficients_list[[as.character(id)]] <- list(Original = original_coefficients, Shifted = shifted_coefficients)
    }
    return(bspline_coefficients_list)
  }
  
  time1<-Sys.time()
  
  # Parameter settings
  seed<-98
  set.seed(seed)
  num_individuals <- 1000
  time_points <- 28
  start_time <- 1
  end_time <- 21
  feature <- "A1"
  alpha<-0.95
  k_range <- 2:8
  timepos_option<-c(4,3,2,1)
  timepos_shift<-c(-1,-2,-3,-4)
  
  data_new_three<- data.frame(ID = integer(), Feature = character(), Time = numeric(), Value = numeric())
  
  # Generate data
  for (i in 1:num_individuals) {
    
    times <- sort(runif(time_points, min = start_time, max = end_time))
    
    if (i <= 300) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 20 + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else if (i > 300 & i <= 600) {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 17 + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else if (i > 600 & i <= 800) {
      
      base_sequence <- 3*sin(0.9*(times+4))
      individual_sequence <- 16+0.5*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    }else {
      times <- sort(runif(round(time_points*1.3,0), min = start_time, max = end_time))
      base_sequence <- 3*sin(0.9*(times+4))
      individual_sequence <- 16+0.5*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
      times<-times*1.3
      
    }
    
    individual_data <- data.frame(ID = rep(i, length(times)), 
                                  Feature = rep(feature, length(times)), 
                                  Time = times, 
                                  Value = individual_sequence)
    data_new_three <- rbind(data_new_three, individual_data)
  }
  
  data_new_three$GroupCategory <- cut(data_new_three$ID,
                                      breaks = c(-Inf, 300, 600, Inf),
                                      labels = c("Group 1", "Group 2", "Group 3"))
  
  data_new_three1<-data_new_three[data_new_three$Time >= start_time & data_new_three$Time <= 17, ]
  p1<-ggplot(data_new_three1, aes(x = Time, y = Value, group = ID, color = GroupCategory)) +
    geom_line(size = 0.2) +
    geom_point(size = 0.5) +
    labs(title = "Scenario 8 True Trajectories", x = "Time", y = "Value", color = "Group") +
    scale_color_manual(values = c("Group 1" = "red", "Group 2" = "blue", "Group 3" ="gold" ),
                       labels = c("Group 1", "Group 2", "Group 3")) +
    theme(
      legend.position = c(0.01, 0.99), 
      legend.justification = c(0.01, 0.99),  
      legend.key.size = unit(0.5, "lines"), 
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 8) 
    )
  
  
  
  
  # Shift time
  data_new_three_shifted <- data_new_three
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(shift_amounts[1:120],rep(0,180),shift_amounts[121:240],rep(0,180),shift_amounts[241:320],rep(0,120),shift_amounts[321:400],rep(0,120) )
  for (i in 1:num_individuals) {
    individual_indices <- data_new_three_shifted$ID == i
    data_new_three_shifted$Time[individual_indices] <- data_new_three_shifted$Time[individual_indices] + shift_amounts[i]
    data_new_three_shifted <- data_new_three_shifted[data_new_three_shifted$Time >= start_time & data_new_three_shifted$Time <= 17 , ]
  }
  
  
  p2<-ggplot(data_new_three_shifted, aes(x = Time, y = Value, group = ID)) +
    geom_line(size = 0.2, color = "skyblue") +
    geom_point(size = 0.5, color = "skyblue") +
    labs(title = "Scenario 8 Observed Trajectories", x = "Time", y = "Value")
  
  source("./EHR-registration-subtype-main/algorithm.R")
  data_recover<-SubtypeAware_Registration(data=data_new_three_shifted,alpha=0.95,k_range=2:8,timepos_option=c(4,3,2,1),tau=0.45,tmin=1,tmax=17,nots=c(8,13))
  
 
  data_recover <- data_recover %>%
    mutate(GroupCategory = factor(GroupCategory, levels = c("Group 1", "Group 2", "Group 3"))) %>%
    arrange(GroupCategory)
  
  p3<-ggplot(data_recover, aes(x = Time, y = Value, group = ID, color = GroupCategory)) +
    geom_line(size = 0.2) +
    geom_point(size = 0.5) +
    labs(title = "Scenario 8 Registered Trajectories", x = "Time", y = "Value", color = "Group") +
    scale_color_manual(values = c("green", "orange","purple")) +
    theme(
      legend.position = c(0.01, 0.99), 
      legend.justification = c(0.01, 0.99),  
      legend.key.size = unit(0.5, "lines"), 
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 8) 
    )

  library(patchwork)
  combined_plot_S8 <- p2 |p3 | p1  
  
  ggsave(filename = "./EHR-registration-subtype-main/Simulation/Figs/NEWFIGS/NEWFig_S8.png",
         plot = combined_plot_S8,
         width = 13,  
         height = 4,  
         dpi = 300    
  )