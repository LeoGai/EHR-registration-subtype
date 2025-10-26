library(readr)
library(tidyverse)
library(cluster)
library(gridExtra)
library(splines)
library(glmnet)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
set.seed(6666)
timepos_option<-c(4,3,2,1)
AKI_reg_base_all_it_10 <- read_csv("./EHR-registration-subtype-main/Real_data/AKI_reg_base_all_it_10.csv")
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

ss3 <- AKI_reg_base_all_it_10 %>%
  group_by(ID) %>%
  slice(1) %>%    
  mutate(TimeShift = Time_reg - Time) %>%
  select(ID, TimeShift) %>%
  ungroup()
id_order <- AKI_data_raw_filtered %>%
  distinct(ID) %>%
  pull(ID)

ss3_sorted <- ss3 %>%
  mutate(ID = factor(ID, levels = id_order)) %>%
  arrange(ID) %>%
  mutate(ID = as.integer(as.character(ID))) 

AKI_data_recover3  <- AKI_data_raw_filtered
shift_recover <- as.numeric(ss3_sorted$TimeShift)

for (i in 1:length(shift_recover)) {
  
  individual_indices <-AKI_data_recover3$ID == unique(AKI_data_recover3$ID)[i]
  AKI_data_recover3$Time[individual_indices] <- AKI_data_recover3$Time[individual_indices] + shift_recover[i]
  AKI_data_recover3 <-AKI_data_recover3[AKI_data_recover3$Time >= 0 & AKI_data_recover3$Time <= 21, ]
}


#######Under 21 days##################
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
AKI_data_recover_21day<-AKI_data_recover3
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
  print("==============k-means recover======================")
  print(round(mean(silhouette_res_recover_21day[, "sil_width"]),3))
  
}

for (j in 2:8) {
  km_res_recover_21day <- pam(data_clustering_recover_21day, k = j)
  silhouette_res_recover_21day <- silhouette(km_res_recover_21day$clustering, dist(data_clustering_recover_21day))
  
  print("=============k-medoids recover=======================")
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



AKI_recover_SC<-AKI_data_recover3[,c(1:4,6)]
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


AKI_data_recover_Clinic<-AKI_data_recover3 %>%
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

