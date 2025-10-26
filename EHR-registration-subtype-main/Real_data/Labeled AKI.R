library(tidyverse)
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


AKI_lab_SC <- AKI_lab_SC %>%
  filter(RECORD_MIN <= 21*24*60)

AKI_lab_SC <- AKI_lab_SC %>%
  group_by(ICU_ID) %>%
  mutate(count = n()) %>%
  filter(count >= 15) %>%
  ungroup()

AKI_lab_SC$BASEVALUE <- AKI_base$Base_SCr

AKI_lab_SC <-AKI_lab_SC %>%
  group_by(ICU_ID) %>%
  filter(
    sum(RECORD_MIN >= 0 & RECORD_MIN <= 7*24*60) >= 3 &   
      sum(RECORD_MIN > 7*24*60 & RECORD_MIN <= 14*24*60) >= 3 &   
      sum(RECORD_MIN > 14*24*60 & RECORD_MIN <= 21*24*60) >= 3  
  ) %>%
  ungroup()

AKI_lab_SC<-AKI_lab_SC[!is.na(AKI_lab_SC[, 4]), ]

AKI_lab_SC_adj <- AKI_lab_SC %>% 
  arrange(ICU_ID, RECORD_MIN) %>%  
  group_by(ICU_ID) %>%  
  filter(RECORD_MIN <= 60*24*7) %>%
  mutate(BASELINE_VALUE = BASEVALUE,  
         BASELINE_TIME = 0, 
         VALUE_RATIO = VALUE / BASELINE_VALUE) %>%  
  mutate(
    MAX_VALUE_2d = ifelse(RECORD_MIN <= 60*24*2, max(VALUE[RECORD_MIN <= 60*24*2], na.rm = TRUE), NA),
    Increase_VALUE_2d=ifelse(RECORD_MIN <= 60*24*2, (max(VALUE[RECORD_MIN <= 60*24*2], na.rm = TRUE)-BASELINE_VALUE), NA),
    MAX_VALUE_RATIO_7d = ifelse(RECORD_MIN <= 60*24*7, max(VALUE_RATIO, na.rm = TRUE), NA)) %>%  
  mutate(SC_CATEGORY = case_when( 
    BASELINE_TIME <= 60*24*2 & (MAX_VALUE_2d >= 4.0) ~ "S4",
    BASELINE_TIME <= 60*24*2 & (MAX_VALUE_RATIO_7d >= 3) ~ "S3",
    BASELINE_TIME <= 60*24*2 & MAX_VALUE_RATIO_7d >= 2 & MAX_VALUE_RATIO_7d < 3 ~ "S2", 
    BASELINE_TIME <= 60*24*2 & MAX_VALUE_RATIO_7d >= 1.5 & MAX_VALUE_RATIO_7d < 2 ~ "S1",
    BASELINE_TIME <= 60*24*2 & (Increase_VALUE_2d >= 0.3) ~ "S0",
    TRUE ~ "unknown" 
  )) %>%
  summarise(SC_CATEGORY = first(SC_CATEGORY)) 

ICU_ID_SC_CATEGORY <- as.matrix(AKI_lab_SC_adj)
ICU_ID_SC_CATEGORY_1118 <- as.data.frame(ICU_ID_SC_CATEGORY)
names(ICU_ID_SC_CATEGORY_1118) <- c("ICU_ID", "SC_CATEGORY")
ICU_ID_SC_CATEGORY_1118$ICU_ID <- as.integer(ICU_ID_SC_CATEGORY_1118$ICU_ID)
saveRDS(ICU_ID_SC_CATEGORY_1118,"./EHR-registration-subtype-main/Real_data/ICU_ID_SC_CATEGORY_1118.rds")
write.csv(ICU_ID_SC_CATEGORY_1118, file = "./EHR-registration-subtype-main/Real_data/ICU_ID_SC_CATEGORY_1118.csv", row.names = FALSE)

