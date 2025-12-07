library(tidyverse)
library(cluster)
library(gridExtra)
library(splines)
library(glmnet)
library(conflicted)
library(patchwork)
library(stringr)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
alpha<-0.95
k_range <- 2:8
timepos_option<-c(4,3,2,1)
set.seed(1111)

AKI_diag<-read.csv("./EHR-registration-subtype-main/Real_data/diagnose_AKI_May_iv.csv")
uni_id<-unique(AKI_diag$ICU_ID)
AKI_lab<-read.csv("./EHR-registration-subtype-main/Real_data/lab_AKI_May_iv.csv")
AKI_base<-read.csv("./EHR-registration-subtype-main/Real_data/AKI_raw_base_0430.csv")
endstage_index <- grepl("End stage renal disease", AKI_diag$DIAGNOSE, ignore.case = TRUE)
endstage_output <- AKI_diag$ICU_ID[endstage_index ]
drop_endstage_id<-unique(endstage_output)
uni_id_01 <- uni_id[!uni_id %in% drop_endstage_id]
burn_index <- grepl("(?i)\\b(?!heartburn)\\w*burn\\b", AKI_diag$DIAGNOSE, perl = TRUE)
burn_output <- AKI_diag$ICU_ID[burn_index ]
uni_id_02 <-uni_id_01 [!uni_id_01  %in% unique(burn_output)]
rd_index <- grepl("renal dialysis", AKI_diag$DIAGNOSE, ignore.case = TRUE)
rd_output <- AKI_diag$ICU_ID[rd_index ]
uni_id_03 <-uni_id_02[!uni_id_02  %in% unique(rd_output)]
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


##################functions###################
compute_max_late <- function(dat) {
  dat %>%
    filter(Time_Bin > 15) %>%            
    group_by(ID, Cluster) %>%
    summarise(
      max_value = max(Value, na.rm = TRUE),  
      .groups = "drop"
    )
}
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

####################recover data####################

AKI_data_recover <- readRDS("./EHR-registration-subtype-main/Real_data/AKI_DA_RE/AKI_data_recover.rds")
timepos_option<-c(1)
j<-6
AKI_data_recover_21day<-AKI_data_recover
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
  coefficients_df0_recover_21day <- rbind(coefficients_df0_recover_21day, data.frame(ID = as.integer(names(bspline_coefficients_list_recover_21day[i])), Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
}

data_clustering_recover_21day <- coefficients_df0_recover_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
rownames(data_clustering_recover_21day)<-coefficients_df0_recover_21day$ID

km_res_recover_21day <- pam(data_clustering_recover_21day, k = j)


AKI_data_recover <- AKI_data_recover |>
  dplyr::mutate(Cluster = km_res_recover_21day[["clustering"]][as.character(ID)])

AKI_data_recover <- AKI_data_recover |>
  dplyr::mutate(
    Cluster = dplyr::recode(
      Cluster,
      `1` = "S[1]^rec",
      `2` = "S[2]^rec",
      `3` = "S[3]^rec",
      `4` = "S[4]^rec",
      `5` = "S[5]^rec",
      `6` = "S[6]^rec"
    )
  )

time_interval <- 1 
AKI_data_recover$Time_Bin <- with(AKI_data_recover, cut(Time, breaks = seq(0,21, by = time_interval), labels = FALSE, include.lowest = TRUE))

AKI_lab_SC_summary <-  AKI_data_recover %>% 
  dplyr::group_by(Cluster, Time_Bin) %>%
  dplyr::summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SD_Value   = sd(Value, na.rm = TRUE),                
    Lower_SD   = Mean_Value - 0.5*SD_Value,
    Upper_SD   = Mean_Value + 0.5*SD_Value,
    .groups = "drop"
  )

draw_key_band_line_point <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(
        col  = NA,
        fill = scales::alpha(data$colour, 0.2)
      )
    ),
    grid::linesGrob(
      x = c(0.1, 0.9), y = 0.5,
      gp = grid::gpar(
        col = data$colour,
        lwd = 0.5
      )
    ),
    grid::pointsGrob(
      x = 0.5, y = 0.5,
      pch = 16,
      size = grid::unit(2, "mm"),
      gp = grid::gpar(
        col = data$colour
      )
    )
  )
}

AKI_lab_SC_summary_early <-AKI_lab_SC_summary%>% filter(Time_Bin <= 10 )
p_recover_early <- ggplot(AKI_lab_SC_summary_early, aes(x = Time_Bin, y = Mean_Value, colour = Cluster, group = Cluster) ) + 
  geom_ribbon( aes(ymin = Lower_SD, ymax = Upper_SD, fill = Cluster), alpha = 0.2, colour = NA, show.legend = FALSE ) + 
  geom_line(size = 0.3, key_glyph = draw_key_band_line_point) + geom_point(size = 0.5, show.legend = FALSE) + 
  labs( title = "Registered: Early SCr Trajectories", x = "Time (days)", y = "SCr Value", colour = "Subtype" ) + 
  scale_colour_discrete( labels = expression( S[1]^"rec", S[2]^"rec", S[3]^"rec", S[4]^"rec", S[5]^"rec", S[6]^"rec" ) ) + 
  theme_minimal(base_size = 14)+ scale_x_continuous(breaks = 1:10)


AKI_lab_SC_summary_mid <-AKI_lab_SC_summary%>% filter(Time_Bin > 10 & Time_Bin <= 15 )

p_recover_mid <- ggplot( AKI_lab_SC_summary_mid, aes(x = Time_Bin, y = Mean_Value, colour = Cluster, group = Cluster) ) + 
  geom_ribbon( aes(ymin = Lower_SD, ymax = Upper_SD, fill = Cluster), alpha = 0.2, colour = NA, show.legend = FALSE ) + 
  geom_line(size = 0.3, key_glyph = draw_key_band_line_point) + geom_point(size = 0.5, show.legend = FALSE) + 
  labs( title = "Registered: Mid SCr Trajectories", x = "Time (days)", y = "SCr Value", colour = "Subtype" ) + 
  scale_colour_discrete( labels = expression( S[1]^"rec", S[2]^"rec", S[3]^"rec", S[4]^"rec", S[5]^"rec", S[6]^"rec" ) ) + 
  theme_minimal(base_size = 14)+ scale_x_continuous(breaks = 11:15)

AKI_lab_SC_summary_late  <- compute_max_late(AKI_data_recover)

p_recover_late<- ggplot(AKI_lab_SC_summary_late,
                        aes(x = Cluster, y = max_value, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  labs(
    title = "Registered: SCr-based Outcomes",
    x = "Subtype",
    y = "Late Peak SCr Value",
    fill = "Subtype" 
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_fill_discrete(labels = function(x) parse(text = x)) +
  theme_minimal(base_size = 14)
####################raw data####################

AKI_data_raw_filtered_21day<-AKI_data_raw_filtered
bspline_coefficients_list_raw_21day<-map_bspline(AKI_data_raw_filtered_21day)
coefficients_df0_raw_21day <- data.frame(ID = integer(), Coef1 = numeric(), Coef2 = numeric(), Coef3 = numeric(), Coef4 = numeric(), Coef5 = numeric(), Coef6 = numeric())

for (i in 1:length(bspline_coefficients_list_raw_21day)) {
  original_coefficients <- bspline_coefficients_list_raw_21day[[i]]$Original
  coef1 <- original_coefficients["(Intercept)"]
  coef2 <- original_coefficients["bsplines1"]
  coef3 <- original_coefficients["bsplines2"]
  coef4 <- original_coefficients["bsplines3"]
  coef5 <- original_coefficients["bsplines4"]
  coef6 <- original_coefficients["bsplines5"]
  coefficients_df0_raw_21day <- rbind(coefficients_df0_raw_21day, data.frame(ID = as.integer(names(bspline_coefficients_list_raw_21day[i])), Coef1 = coef1, Coef2 = coef2, Coef3 = coef3, Coef4 = coef4, Coef5 = coef5, Coef6 = coef6))
}

data_clustering_raw_21day <- coefficients_df0_raw_21day[, c("Coef1", "Coef2", "Coef3", "Coef4", "Coef5", "Coef6")]
rownames(data_clustering_raw_21day)<-coefficients_df0_raw_21day$ID

km_res_raw_21day <- pam(data_clustering_raw_21day, k = j)


AKI_data_raw_filtered <- AKI_data_raw_filtered |>
  dplyr::mutate(Cluster = km_res_raw_21day[["clustering"]][as.character(ID)])

AKI_data_raw_filtered <- AKI_data_raw_filtered |>
  dplyr::mutate(
    Cluster = dplyr::recode(
      Cluster,
      `1` = "S[1]^raw",
      `2` = "S[2]^raw",
      `3` = "S[3]^raw",
      `4` = "S[4]^raw",
      `5` = "S[5]^raw",
      `6` = "S[6]^raw"
    )
  )

time_interval <- 1 
AKI_data_raw_filtered$Time_Bin <- with(AKI_data_raw_filtered, cut(Time, breaks = seq(0,21, by = time_interval), labels = FALSE, include.lowest = TRUE))

AKI_lab_SC_summary_raw <- AKI_data_raw_filtered %>% 
  dplyr::group_by(Cluster, Time_Bin) %>%
  dplyr::summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SD_Value   = sd(Value, na.rm = TRUE),                 
    Lower_SD   = Mean_Value - 0.5*SD_Value,
    Upper_SD   = Mean_Value + 0.5*SD_Value,
    .groups = "drop"
  )



AKI_lab_SC_summary_early_raw <-AKI_lab_SC_summary_raw %>% filter(Time_Bin <= 10 )
p_recover_early_raw <- ggplot(AKI_lab_SC_summary_early_raw, aes(x = Time_Bin, y = Mean_Value, colour = Cluster, group = Cluster) ) + 
  geom_ribbon( aes(ymin = Lower_SD, ymax = Upper_SD, fill = Cluster), alpha = 0.2, colour = NA, show.legend = FALSE ) + 
  geom_line(size = 0.3, key_glyph = draw_key_band_line_point) + geom_point(size = 0.5, show.legend = FALSE) + 
  labs( title = "Raw: Early SCr Trajectories", x = "Time (days)", y = "SCr Value", colour = "Subtype" ) + 
  scale_colour_discrete( labels = expression( S[1]^"raw", S[2]^"raw", S[3]^"raw", S[4]^"raw", S[5]^"raw", S[6]^"raw" ) ) + 
  theme_minimal(base_size = 14)+ scale_x_continuous(breaks = 1:10)


AKI_lab_SC_summary_mid_raw <-AKI_lab_SC_summary_raw %>% filter(Time_Bin > 10 & Time_Bin <= 15 )

p_recover_mid_raw <- ggplot( AKI_lab_SC_summary_mid_raw, aes(x = Time_Bin, y = Mean_Value, colour = Cluster, group = Cluster) ) + 
  geom_ribbon( aes(ymin = Lower_SD, ymax = Upper_SD, fill = Cluster), alpha = 0.2, colour = NA, show.legend = FALSE ) + 
  geom_line(size = 0.3, key_glyph = draw_key_band_line_point) + geom_point(size = 0.5, show.legend = FALSE) + 
  labs( title = "Raw: Mid SCr Trajectories", x = "Time (days)", y = "SCr Value", colour = "Subtype" ) + 
  scale_colour_discrete( labels = expression( S[1]^"raw", S[2]^"raw", S[3]^"raw", S[4]^"raw", S[5]^"raw", S[6]^"raw" ) ) + 
  theme_minimal(base_size = 14)+ scale_x_continuous(breaks = 11:15)

AKI_lab_SC_summary_late_raw  <- compute_max_late(AKI_data_raw_filtered)

p_recover_late_raw<- ggplot(AKI_lab_SC_summary_late_raw,
                        aes(x = Cluster, y = max_value, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  labs(
    title = "Raw: SCr-based Outcomes",
    x = "Subtype",
    y = "Late Peak SCr Value",
    fill = "Subtype" 
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_fill_discrete(labels = function(x) parse(text = x)) +
  theme_minimal(base_size = 14)

combined_plot_subda <- 
  (p_recover_early      | p_recover_mid      | p_recover_late) /
  (p_recover_early_raw  | p_recover_mid_raw  | p_recover_late_raw)


ggsave(filename = "./EHR-registration-subtype-main/Real_data/AKI_DA_RE/combined_plot_subda.png",
       plot = combined_plot_subda,
       width = 16,  
       height = 9,  
       dpi = 300    
)

###############Diagnoses##########################################

pick_col <- function(df, cand) {
  nms <- names(df); hit <- NULL
  for (x in cand) {
    idx <- which(tolower(nms) == tolower(x))
    if (length(idx)) { 
      hit <- nms[idx[1]]
      break
    }
  }
  if (is.null(hit)) {
    stop(sprintf(
      "Cannot find any of the candidate columns: %s. Available columns: %s",
      paste(cand, collapse = ", "),
      paste(nms, collapse = ", ")
    ))
  }
  hit
}

norm_txt <- function(x) x |>
  stringr::str_to_lower() |>
  stringr::str_replace_all("[^a-z0-9+ ]"," ") |>
  stringr::str_squish()

make_patterns <- function(include_urosepsis = TRUE, include_bacteremia = FALSE){
  
  # drop AKI/ARF
  pat_excl_aki <- regex("(acute\\s+(kidney|renal)\\s+(failure|injury))|\\baki\\b", TRUE)
  
  # 1) Hyperkalemia
  pat_hyperk <- regex(
    "hyperkal|hyperpotass|\\bk\\s*\\+\\s*(high|elevat)|high potassium|elevated potassium|\\bhyper\\s?k\\b",
    TRUE
  )
  neg_hyperk <- regex("hypokal|low potassium|pseudo\\s*hyperkal", TRUE)
  
  # 2) Hyponatremia / Hyposmolality
  pat_hypona <- regex("hyponatrem|hyposmolal", TRUE)
  
  # 3) Acute respiratory failure
  pat_resp <- regex("\\bacute\\s+respiratory\\s+failure(\\s+with\\s+hypoxia)?\\b", TRUE)
  
  # 4) Sepsis spectrum
  pat_sepsis <- regex(paste0(
    "septic shock|severe sepsis|septic(a)?emia|\\bsepsis\\b",
    if (include_urosepsis)  "|\\burosepsis\\b" else "",
    if (include_bacteremia) "|\\bbacteremia\\b" else ""
  ), TRUE)
  neg_sepsis <- regex("asepsis|antisepsis|screening.*sepsis|history of sepsis|rule out sepsis|sepsis risk", TRUE)
  
  # 5) Essential hypertension
  pat_htn <- regex(
    "essential\\s+.*hypertension|unspecified essential hypertension",
    TRUE
  )
  
  # 6) Hyperlipidemia
  pat_hld <- regex("hyperlipidemia", TRUE)
  
  # 7) Atrial fibrillation
  pat_afib <- regex("atrial fibrillation", TRUE)
  
  # 8) Ascites
  pat_ascites <- regex("\\bascites\\b", TRUE)
  
  # 9) UTI
  pat_uti <- regex("urinary tract infection|\\buti\\b", TRUE)
  
  # 10) CKD
  pat_ckd <- regex("chronic kidney disease|\\bckd\\b|chronic renal (insufficiency|failure)", TRUE)
  
  list(
    pat_excl_aki = pat_excl_aki,
    pat_hyperk   = pat_hyperk,
    neg_hyperk   = neg_hyperk,
    pat_hypona   = pat_hypona,
    pat_resp     = pat_resp,
    pat_sepsis   = pat_sepsis,
    neg_sepsis   = neg_sepsis,
    pat_htn      = pat_htn,
    pat_hld      = pat_hld,
    pat_afib     = pat_afib,
    pat_ascites  = pat_ascites,
    pat_uti      = pat_uti,
    pat_ckd      = pat_ckd
  )
}

merge_count_plot <- function(AKI_data_recover, AKI_diag,
                             cluster_name = "Subtype 5",
                             top_n = 20,
                             include_urosepsis = TRUE,
                             include_bacteremia = FALSE){
  
  cl_col  <- pick_col(AKI_data_recover, c("Cluster"))
  id_rec  <- pick_col(AKI_data_recover, c("ID"))
  id_diag <- pick_col(AKI_diag,         c("ICU_ID"))
  dx_col  <- pick_col(AKI_diag,         c("DIAGNOSE"))
  
  ids <- AKI_data_recover %>%
    filter(.data[[cl_col]] == cluster_name) %>%
    pull(all_of(id_rec)) %>% unique()
  N <- length(ids)

  diag_df <- AKI_diag %>%
    filter(.data[[id_diag]] %in% ids, !is.na(.data[[dx_col]])) %>%
    transmute(ICU_ID = .data[[id_diag]],
              DIAGNOSE = .data[[dx_col]],
              dx_norm  = norm_txt(.data[[dx_col]]))
  
  P <- make_patterns(include_urosepsis, include_bacteremia)
  
  diag_merge <- diag_df %>%
    filter(!str_detect(dx_norm, P$pat_excl_aki)) %>%
    mutate(
      merged_concept = case_when(
        str_detect(dx_norm, P$pat_hyperk) & !str_detect(dx_norm, P$neg_hyperk) ~ "Hyperkalemia (merged)",
        str_detect(dx_norm, P$pat_hypona)                                      ~ "Hyponatremia / Hyposmolality (merged)",
        str_detect(dx_norm, P$pat_resp)                                        ~ "Acute respiratory failure (merged)",
        str_detect(dx_norm, P$pat_sepsis) & !str_detect(dx_norm, P$neg_sepsis) ~ "Sepsis spectrum (merged)",
        str_detect(dx_norm, P$pat_htn)     ~ "Essential hypertension (merged)",
        str_detect(dx_norm, P$pat_hld)     ~ "Hyperlipidemia (merged)",
        str_detect(dx_norm, P$pat_afib)    ~ "Atrial fibrillation (merged)",
        str_detect(dx_norm, P$pat_ascites) ~ "Ascites (merged)",
        str_detect(dx_norm, P$pat_uti)     ~ "Urinary tract infection (merged)",
        str_detect(dx_norm, P$pat_ckd)     ~ "Chronic kidney disease (merged)",
        
        TRUE ~ NA_character_
      ),
      label_merged = if_else(is.na(merged_concept), DIAGNOSE, merged_concept)
    ) %>%
    select(ICU_ID, label_merged)
  
  freq_tbl <- diag_merge %>%
    distinct(ICU_ID, label_merged) %>%
    count(label_merged, name = "n_patients", sort = TRUE) %>%
    mutate(per100 = 100 * n_patients / N)
  
  # Top N
  plot_df <- freq_tbl %>%
    slice_max(n_patients, n = top_n) %>%
    mutate(label_show = stringr::str_wrap(label_merged, width = 40))
  
  cluster_expr <- tryCatch(
    parse(text = cluster_name)[[1]],
    error = function(e) cluster_name
  )
  
  p <- ggplot(plot_df,
              aes(x = per100, y = reorder(label_show, per100))) +
    geom_col(width = 0.7, fill = "#00B2EE") +   
    geom_text(
      aes(label = sprintf("%.1f%%", per100)),  
      hjust = -0.15,
      size  = 3.2,
      color = "#333333"
    ) +
    scale_x_continuous(
      labels = function(x) sprintf("%.0f%%", x),          
      expand = expansion(mult = c(0, 0.25))
    ) +
    labs(
      x = bquote("Patients with diagnosis (% of " * .(cluster_expr) * " patients)"),
      y = NULL,
      title = bquote(.(cluster_expr) * " cohort: Top " * .(top_n) *
                       " diagnoses by patient proportion")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 9),
      plot.title         = element_text(face = "bold", hjust = 0),
      axis.title.x       = element_text(margin = margin(t = 6)),
      plot.margin        = margin(t = 5, r = 15, b = 5, l = 5)
    )
  
  list(table = freq_tbl, plot = p, N = N)
}




res_S3 <- merge_count_plot(AKI_data_recover, AKI_diag, "S[3]^rec", top_n = 20)
#res_S3$plot


res_S5 <- merge_count_plot(AKI_data_recover, AKI_diag, "S[5]^rec", top_n = 20)
#res_S5$plot

res_S6 <- merge_count_plot(AKI_data_recover, AKI_diag, "S[6]^rec", top_n = 20)
#res_S6$plot



p_S3 <- res_S3$plot +
  labs(
    title = expression(S[3]^rec),
    x = "Patients with diagnosis (%)",
    y = "Diagnosis"
  )

p_S5 <- res_S5$plot +
  labs(
    title = expression(S[5]^rec),
    x = "Patients with diagnosis (%)",
    y = NULL
  )

p_S6 <- res_S6$plot +
  labs(
    title = expression(S[6]^rec),
    x = "Patients with diagnosis (%)",
    y = NULL
  )

combined_diag <- p_S3 | p_S5 | p_S6

combined_diag <- combined_diag +
  plot_annotation(
    title = "Registered: Top diagnoses by patient proportion for selected subtypes"
  )

combined_diag <- combined_diag +
  ggplot2::theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  )


ggsave(filename = "./EHR-registration-subtype-main/Real_data/AKI_DA_RE/combined_diag.png",
       plot = combined_diag ,
       width = 16,  
       height = 6,  
       dpi = 300    
)


