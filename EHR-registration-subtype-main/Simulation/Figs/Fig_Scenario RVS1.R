rm(list = ls())
library(tidyverse)    
library(cluster)    
library(splines)
library(glmnet)
library(ggplot2)
library(gridExtra)
library(patchwork)
  
  
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
  data_new_four1<-data_new_four1[data_new_four1$ID %in% c(7,89,265,353,575,617,837,900), ]
  
  
  p1<-ggplot(data_new_four1, aes(x = Time, y = Value, group = ID, color = GroupCategory)) +
    geom_line(size = 0.2) +
    geom_point(size = 1.5) +
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
 
  data_new_four_shifted1<-data_new_four_shifted[data_new_four_shifted$ID %in% c(7,89,265,353,575,617,837,900), ]
  
  
  
  
  
  p2<-ggplot(data_new_four_shifted1, aes(x = Time, y = Value, group = ID)) +
    geom_line(size = 0.2, color = "skyblue") +
    geom_point(size = 1.5, color = "skyblue") +
    labs(title = "Scenario 4 Observed Trajectories", x = "Time", y = "Value")
  
  
  
  
 
library(patchwork)
combined_plot_RVS1 <- p2 | p1

ggsave(filename = "./EHR-registration-subtype-main/Simulation/Figs/NEWFIGS/NEWFig_RVS1.png",
       plot = combined_plot_RVS1,
       width = 10,  
       height = 4,  
       dpi = 300    
)
