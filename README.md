# Experiments for Subtype-Aware Registration of Longitudinal Electronic Health Records

This repository contains the code scripts used in the following paper:

- **Title:** Subtype-Aware Registration of Longitudinal Electronic Health Records
- **Authors:** Xin Gai, Shiyi Jiang, and Anru R. Zhang

If you have any question about this code, please contact Xin Gai (xg.biostat@gmail.com).

# Reproducing Experiments

## Setup
The experiments require:

- **R version**: 4.5.1
- **R libraries**:

  - `doParallel` : 1.0.17  
  - `parallel` : 4.5.1  
  - `cluster` : 2.1.8.1  
  - `ggplot2` : 4.0.0  
  - `gridExtra` : 2.3
  - `patchwork` : 1.3.2  
  - `splines` : 4.5.1  
  - `glmnet` : 4.1.10  
  - `tidyverse` : 2.0.0  
  - `kernlab` : 0.9.33  
  - `dtwclust` : 6.0.0  
  - `proxy` : 0.4.27  
  - `mclust` : 6.1.1  
  - `aricode` : 1.0.3  
  - `clue` : 0.3.66 
  
- **Python version**: 3.6.8
- **Python libraries**:

  - `pandas`: 1.1.5
  - `numpy`: 1.18.0
  - `scipy`: 1.4.1
  - `scikit-learn`: 0.24.2
  - `torch`: 1.4.0
  - `matplotlib`: 3.3.3
  - `statsmodels`: 0.12.2

## Contents

* `algorithm.R` contains the main algorithm discussed in the paper.

### Simulation Study

- The `Simulation/` directory contains all files related to simulation studies. This directory includes eight R scripts: `Scenario1.R`, `Scenario2.R`, `Scenario3.R`, `Scenario4.R`, `Scenario5.R`, `Scenario6.R`, `Scenario7.R`, and `Scenario8.R`. Each script corresponds to a specific simulation scenario (Scenario 1-8) used in the study. Running these scripts will generate the results summarized in **Table 1** of the manuscript.

- The `Simulation/Figs/` directory contains all R scripts for visualizing the simulation study results. Running these scripts will generate the eight Simulation Results Figures (Figures 1 and 9–15) included in the manuscript.

#### Ablation Studies
- The `Simulation_ablation1/`, `Simulation_ablation2/`, and `Simulation_ablation3/` directories hold scripts for the smoothing-spline (**Table 7**), RKHS (**Table 8**), and k-means (**Table 9**) ablation studies, respectively.

#### Sensitivity Analyses
- The `Sensitivity_analyses/` folder holds six sub-folders: `SA_1_3/`, `SA_4_6/`, `SA_7_9/`, `SA_10_14/`, `SA_15_21/`, and `SA_22_25/`. 
- Each sub-folder contains a script of the same name (e.g., SA_1_3.R) that conducts sensitivity analysis on one specific aspect: outliers, missing data, noise, $\alpha$, M and $\tau$. After running the individual scripts, execute SA_Figs.R to create Figures 3–8 for the manuscript.

#### Additional Comparison Results
- The `Comparison_DTW/`directory contains scripts for comparing our method with the DTW method on the clustering task (**Table 10**).
  
### Real Data Experiment

- The `Real_data/` directory contains all files on real data analysis.

#### Key Files and Directories
- **`Labeled_AKI.R`**: Generates AKI classification labels, saving them as `ICU_ID_SC_CATEGORY_1118.rds`.

- **`AKI_DA_RE/`**: Contains scripts for implementing the Registration algorithm and performing related tasks (Unsupervised Learning, Clinical Evaluation) using **k-medoids** clustering.  
  - **`AKI_DA.R`**: Executes the Registration algorithm for K values ranging from [2, 8], generating the **Registered Data** saved as `AKI_data_recover.csv`.  
  - **`AKI_DA_LR.R`**: Executes the Registration algorithm for K values ranging from [5, 8].  
  - **`AKI_DA_SR.R`**: Executes the Registration algorithm for K values ranging from [2, 4].
  - **`Phenotypefigs.R`**: Generates the two figures in Section 4.2.1.  
  Running these scripts generates the results corresponding to **Tables 2, 3, and 4** in the manuscript.

- **`AKI_DA_KM_RE/`**: Contains scripts for implementing the Registration algorithm and performing related tasks (Unsupervised Learning, Clinical Evaluation) using **k-means** clustering.  
  - **`AKI_DA_KM.R`**: Executes the Registration algorithm for K values ranging from [2, 8], generating the **Registered Data** saved as `AKI_data_recover2.csv`.  
  - **`AKI_DA_KM_LR.R`**: Executes the Registration algorithm for K values ranging from [5, 8].  
  - **`AKI_DA_KM_SR.R`**: Executes the Registration algorithm for K values ranging from [2, 4].  
  Running these scripts generates the results corresponding to **Tables 2, 3, and 4** in the manuscript.

- **`real_data_prediction/mimic_iv_data_pull.py`**: Extracts all raw Creatinine records of the desired cohort from MIMIC-IV dataset in the format of csv.
- **`real_data_prediction/gen_AKI_labels.py`**: Generates baseline Creatinine levels and corresponding severe AKI labels for subjects.
- **`real_data_prediction/AKI_predict_main.py`**: Reports severe AKI prediction metrics (including metrics from subgroup analysis on age and gender) in terms of mean and 95% confidence interval (**Tables 5, 6, 7, and 8**).
- **`real_data_prediction/predict_utils.py`**: Contains helper functions for performing severe AKI prediction and subgroup analysis on age and gender.
- **`real_data_prediction/registration_baseline/gen_reg_data.py`**: Generates recovered data using the baseline method.

#### Data Source
The dataset used for real data analysis is derived from the **MIMIC-IV Database**. For more details, visit the [MIMIC-IV database](https://physionet.org/content/mimiciv/3.1/).
