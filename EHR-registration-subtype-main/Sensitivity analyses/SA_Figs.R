rm(list = ls())
library(patchwork)
library(ggplot2)
library(tidyverse) 

Reg_SA1 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_1_3/Reg_SA1.rds")
Reg_SA2 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_1_3/Reg_SA2.rds")
Reg_SA3 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_1_3/Reg_SA3.rds")
Reg_SA4 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_4_6/Reg_SA4.rds")
Reg_SA5 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_4_6/Reg_SA5.rds")
Reg_SA6 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_4_6/Reg_SA6.rds")
Reg_SA7 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_7_9/Reg_SA7.rds")
Reg_SA8 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_7_9/Reg_SA8.rds")
Reg_SA9 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_7_9/Reg_SA9.rds")
Reg_SA0 <- readRDS("./EHR-registration-subtype-main/Simulation/FReg_Simu3.rds")


#### 1. Outliers Scenario ####

df_outlier <- rbind(
  cbind(Scenario = "SA0", as.data.frame(Reg_SA0)),
  cbind(Scenario = "100", as.data.frame(Reg_SA1)),
  cbind(Scenario = "200", as.data.frame(Reg_SA2)),
  cbind(Scenario = "300", as.data.frame(Reg_SA3))
)
df_outlier_long <- df_outlier %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_outlier_long$Scenario <- factor(
  df_outlier_long$Scenario,
  levels=c("SA0","100","200","300"),
  ordered=TRUE
)
summary_data <- df_outlier_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=2,shape=21,fill="white") +
    scale_x_continuous(breaks=unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x = "Number of Outliers", y = "Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}


p_outlier1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0,1)) +
  ggtitle(metric_labels["ACC"])
p_outlier2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_outlier3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0,1)) +
  ggtitle(metric_labels["MAE"])

p_outlier <- p_outlier1| p_outlier2 | p_outlier3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", outlier)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"   
    )
  )

print(p_outlier)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_outlier.png",
       plot = p_outlier,
       width = 10,  
       height = 6,  
       dpi = 500    
)

#### 2. Missing Data Scenario  ####

df_missing <- rbind(
  cbind(Scenario = "SA0", as.data.frame(Reg_SA0)),
  cbind(Scenario = "10%", as.data.frame(Reg_SA4)),
  cbind(Scenario = "15%", as.data.frame(Reg_SA5)),
  cbind(Scenario = "20%", as.data.frame(Reg_SA6))
)
df_missing_long <- df_missing %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_missing_long$Scenario <- factor(
  df_missing_long$Scenario,
  levels=c("SA0","10%","15%","20%"),
  ordered=TRUE
)
summary_data <- df_missing_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=2,shape=21,fill="white") +
    scale_x_continuous(breaks = unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x = "Missing Rate", y = "Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

p_missing1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0.6,1)) +
  ggtitle(metric_labels["ACC"])
p_missing2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0.8,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_missing3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0,1)) +
  ggtitle(metric_labels["MAE"])

p_missing <- p_missing1| p_missing2 | p_missing3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", missing)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"  
    )
  )

print(p_missing)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_missing.png",
       plot = p_missing,
       width = 10,  
       height = 6,  
       dpi = 500    
)

#### 3. Increased Error Variance Scenario  ####

df_noise <- rbind(
  cbind(Scenario = "SA0", as.data.frame(Reg_SA0)),
  cbind(Scenario = "15%", as.data.frame(Reg_SA7)),
  cbind(Scenario = "30%", as.data.frame(Reg_SA8)),
  cbind(Scenario = "45%", as.data.frame(Reg_SA9))
)
df_noise_long <- df_noise %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_noise_long$Scenario <- factor(
  df_noise_long$Scenario,
  levels=c("SA0","15%","30%","45%"),
  ordered=TRUE
)
summary_data <- df_noise_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=2,shape=21,fill="white") +
    scale_x_continuous(breaks = unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x = "Relative Noise Increase", y = "Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

p_noise1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0.6,1)) +
  ggtitle(metric_labels["ACC"])
p_noise2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0.8,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_noise3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0,1)) +
  ggtitle(metric_labels["MAE"])

p_noise <- p_noise1| p_noise2 | p_noise3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", noise)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"   
    )
  )

print(p_noise)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_noise.png",
       plot = p_noise,
       width = 10,  
       height = 6,  
       dpi = 500    
)




#### 4. alpha Scenario ####

Reg_SA10 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_10_14/FReg_SA10.rds")
Reg_SA11 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_10_14/FReg_SA11.rds")
Reg_SA12 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_10_14/FReg_SA12.rds")
Reg_SA13 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_10_14/FReg_SA13.rds")
Reg_SA14 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_10_14/FReg_SA14.rds")

df_alpha <- rbind(
  cbind(Scenario = "0.86", as.data.frame(Reg_SA14)),
  cbind(Scenario = "0.89", as.data.frame(Reg_SA13)),
  cbind(Scenario = "0.92", as.data.frame(Reg_SA12)),
  cbind(Scenario = "0.95", as.data.frame(Reg_SA11)),
  cbind(Scenario = "0.98", as.data.frame(Reg_SA10))
)
df_alpha_long <- df_alpha %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_alpha_long$Scenario <- factor(
  df_alpha_long$Scenario,
  levels=c("0.86","0.89","0.92","0.95","0.98"),
  ordered=TRUE
)
summary_data <- df_alpha_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=3,shape=21,fill="white") +
    scale_x_continuous(breaks=unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x=expression(alpha), y="Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

p_alpha1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0.4,1)) +
  ggtitle(metric_labels["ACC"])
p_alpha2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0.8,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_alpha3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0.1,1.0)) +
  ggtitle(metric_labels["MAE"])

p_alpha <- p_alpha1| p_alpha2 | p_alpha3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", alpha)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"  
    )
  )

print(p_alpha)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_alpha.png",
       plot = p_alpha,
       width = 10,  
       height = 6,  
       dpi = 500    
)



#### 5. M Scenario ####

Reg_SA15 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA15.rds")
Reg_SA16 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA16.rds")
Reg_SA17 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA17.rds")
Reg_SA18 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA18.rds")
Reg_SA19 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA19.rds")
Reg_SA20 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA20.rds")
Reg_SA21 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_15_21/Reg_SA21.rds")


df_M <- rbind(
  cbind(Scenario = "2", as.data.frame(Reg_SA15)),
  cbind(Scenario = "3", as.data.frame(Reg_SA16)),
  cbind(Scenario = "4", as.data.frame(Reg_SA17)),
  cbind(Scenario = "5", as.data.frame(Reg_SA18)),
  cbind(Scenario = "6", as.data.frame(Reg_SA19)),
  cbind(Scenario = "8", as.data.frame(Reg_SA20)),
  cbind(Scenario = "10", as.data.frame(Reg_SA21))
)

df_M_long <- df_M %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_M_long$Scenario <- factor(
  df_M_long$Scenario,
  levels=c("2","3","4","5","6","8","10"),
  ordered=TRUE
)
summary_data <- df_M_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=3,shape=21,fill="white") +
    scale_x_continuous(breaks=unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x=expression(M), y="Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}


p_M1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0.7,0.8)) +
  ggtitle(metric_labels["ACC"])
p_M2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0.95,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_M3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0.2,0.4)) +
  ggtitle(metric_labels["MAE"])

p_M <- p_M1| p_M2 | p_M3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", M)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"   
    )
  )

print(p_M)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_M.png",
       plot = p_M,
       width = 10,  
       height = 6,  
       dpi = 500    
)


#####tau###################

Reg_SA22 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_22_25/Reg_SA22.rds")
Reg_SA23 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_22_25/Reg_SA23.rds")
Reg_SA24 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_22_25/Reg_SA24.rds")
Reg_SA25 <- readRDS("./EHR-registration-subtype-main/Sensitivity analyses/SA_22_25/Reg_SA25.rds")


df_tau <- rbind(
  cbind(Scenario = "0.45", as.data.frame(Reg_SA22)),
  cbind(Scenario = "0.48", as.data.frame(Reg_SA23)),
  cbind(Scenario = "0.51", as.data.frame(Reg_SA24)),
  cbind(Scenario = "0.54", as.data.frame(Reg_SA25))
)


df_tau_long <- df_tau %>%
  pivot_longer(c("ACC","LQ1","MAE"), names_to="Metric", values_to="Value")
df_tau_long$Scenario <- factor(
  df_tau_long$Scenario,
  levels=c("0.45","0.48","0.51","0.54"),
  ordered=TRUE
)
summary_data <- df_tau_long %>%
  group_by(Scenario, Metric) %>%
  summarize(
    Q2.5 = quantile(Value,0.25,na.rm=TRUE),
    #Mean = mean(Value,na.rm=TRUE),
    Mean =quantile(Value,0.5,na.rm=TRUE),
    Q97.5 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    Scenario_num = as.numeric(Scenario),
    x_min = Scenario_num - 0.15,
    x_max = Scenario_num + 0.15
  )

metric_labels <- c(ACC="Rate of Exact Recovery",
                   LQ1="Rate of Recovery Error ≤ 1",
                   MAE="MAE (days)")
metric_colors <- c(ACC="dodgerblue",LQ1="seagreen",MAE="orange")
unique_breaks <- sort(unique(summary_data$Scenario_num))
unique_labels <- levels(summary_data$Scenario)

make_plot <- function(df, ylim_range){
  ggplot(df, aes(x=Scenario_num)) +
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=Q2.5,ymax=Q97.5,fill=Metric),
              alpha=0.3,color=NA) +
    geom_line(aes(y=Mean,color=Metric),size=1) +
    geom_point(aes(y=Mean,color=Metric),size=2,shape=21,fill="white") +
    scale_x_continuous(breaks = unique_breaks, labels = unique_labels) +
    scale_fill_manual(values=metric_colors) +
    scale_color_manual(values=metric_colors) +
    labs(x=expression(tau), y = "Value") +
    coord_cartesian(ylim=ylim_range) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}


p_tau1 <- make_plot(filter(summary_data, Metric=="ACC"), c(0,1)) +
  ggtitle(metric_labels["ACC"])
p_tau2 <- make_plot(filter(summary_data, Metric=="LQ1"), c(0.9,1.0)) +
  ggtitle(metric_labels["LQ1"])
p_tau3 <- make_plot(filter(summary_data, Metric=="MAE"), c(0,1)) +
  ggtitle(metric_labels["MAE"])

p_tau <- p_tau1| p_tau2 | p_tau3 + 
  plot_annotation(
    title = expression(paste("Sensitivity Analysis: ", tau)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.title.position = "plot"  
    )
  )

print(p_tau)


ggsave(filename = "./EHR-registration-subtype-main/Sensitivity analyses/SA_FIGS/SAFig_tau.png",
       plot = p_tau,
       width = 10,  
       height = 6,  
       dpi = 500    
)