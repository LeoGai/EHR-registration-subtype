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
  seed<-485
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
  
 
  data_new_four <- data.frame(ID = integer(), Feature = character(), Time = numeric(), Value = numeric())
  
  # Generate data
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
      individual_sequence <- 38-0.5*times + base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    } else {
      base_sequence <- 3*sin(0.6*(times+4))
      individual_sequence <- 21 +0.9*times+ base_sequence + rnorm(length(times), mean = 0, sd =0.8)
    }
    
    individual_data <- data.frame(ID = rep(i, length(times)), 
                                  Feature = rep(feature, length(times)), 
                                  Time = times, 
                                  Value = individual_sequence)
    data_new_four <- rbind(data_new_four, individual_data)
  }
  
  
  
  data_new_four$GroupCategory <- cut(data_new_four$ID,
                                     breaks = c(-Inf, 250,500, 750, Inf),
                                     labels = c("Group 1", "Group 2", "Group 3","Group 4"))
  
  
  data_new_four1<-data_new_four[data_new_four$Time >= start_time & data_new_four$Time <= 17, ]
  p1<-ggplot(data_new_four1, aes(x = Time, y = Value, group = ID, color = GroupCategory)) +
    geom_line(size = 0.2) +
    geom_point(size = 0.5) +
    labs(title = "Scenario 4 True Trajectories", x = "Time", y = "Value", color = "Group") +
    scale_color_manual(values = c("Group 1" = "red", "Group 2" = "blue", "Group 3" = "gold", "Group 4" = "orchid"),
                       labels = c("Group 1", "Group 2", "Group 3","Group 4")) +
    theme(
      legend.position = c(0.99, 1), # Position the legend inside the plot
      legend.justification = c(0.99, 1), # Right and top justify the legend
      legend.key.size = unit(0.5, "lines"), # Adjust the size of the legend keys
      legend.text = element_text(size = 8), # Adjust the size of the legend text
      legend.title = element_text(size = 8) # Adjust the size of the legend title
    )
  

  
  # Shift time
  data_new_four_shifted <- data_new_four
  set.seed(seed)
  shift_amounts <- sample(timepos_shift, 400, replace = TRUE)
  shift_amounts<-c(shift_amounts[1:100],rep(0,150),shift_amounts[101:200],rep(0,150),shift_amounts[201:300],rep(0,150),shift_amounts[301:400],rep(0,150) )
  for (i in 1:num_individuals) {
    individual_indices <- data_new_four_shifted$ID == i
    data_new_four_shifted$Time[individual_indices] <- data_new_four_shifted$Time[individual_indices] + shift_amounts[i]
    data_new_four_shifted <- data_new_four_shifted[data_new_four_shifted$Time >= start_time  & data_new_four_shifted$Time <= 17 , ]
  }
 
  p2<-ggplot(data_new_four_shifted, aes(x = Time, y = Value, group = ID)) +
    geom_line(size = 0.2, color = "skyblue") +
    geom_point(size = 0.5, color = "skyblue") +
    labs(title = "Scenario 4 Observed Trajectories", x = "Time", y = "Value")
  
  
  source("./EHR-registration-subtype-main/algorithm.R")
  data_recover<-SubtypeAware_Registration(data=data_new_four_shifted,alpha=0.95,k_range=2:8,timepos_option=c(4,3,2,1),tau=0.45,tmin=1,tmax=17,nots=c(8,13))
  
  
 
  
  p3<-ggplot(data_recover, aes(x = Time, y = Value, group = ID, color = GroupCategory)) +
    geom_line(size = 0.2) +
    geom_point(size = 0.5) +
    labs(title = "Scenario 4 Registered Trajectories", x = "Time", y = "Value", color = "Group")  +
    scale_color_manual(values = c("green", "orange","purple","coral")) +
    theme(
      legend.position = c(0.99, 1), 
      legend.justification = c(0.99, 1), 
      legend.key.size = unit(0.5, "lines"), 
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 8) 
    )
 

combined_plot_S4 <- p2 |p3 | p1

ggsave(filename = "./EHR-registration-subtype-main/Simulation/Figs/NEWFIGS/NEWFig_S4.png",
       plot = combined_plot_S4,
       width = 13,  
       height = 4,  
       dpi = 300    
)
