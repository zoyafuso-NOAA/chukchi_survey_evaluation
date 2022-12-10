## Workflow:

### 1) [Data Inputs](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/data_processing)

* [Process_AK_BTS_data_arctic.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/data_processing/Process_AK_BTS_data_arctic.R): Pulls Chukchi 1990 and 2012 otter (83-112 eastern bottom trawl) catch and effort data from RACE database.
* [Process_AK_IERL_data_arctic.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/data_processing/Process_AK_IERL_data_arctic.R): Pulls Chukchi 2012, 2017, and 2019 beam trawl data from IERL survey.
* [aggregates_species.csv](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/data/fish_data/aggregates_species.csv): table defining species taxon groups, called in the two above scripts. 

### 2) [Operating Model](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/VAST_fitting)

* [VAST_fitting_general.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/VAST_fitting/VAST_fitting_general.R): Fits single-species VAST models with different model configurations (e.g., spatial/spatiotemporal fields, observation models). Saves the best (lowest-aic) model, runs diagnostics, and simulates densities.  
* [optimization_data.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/survey_optimization/optimization_data.R): Sets up constants used in the (subsequent) STRS design optimization and survey simulation scripts.

### 3) [Stratified Survey Optimization](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/survey_optimization)
* [combined_ss_ms_optimization.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/survey_optimization/combined_ss_ms_optimization.R): Conducts stratified-random design optimization for both gears.

### 4) [Survey Simulation](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/survey_optimization/) 
* [simulate_survey.R](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/blob/main/code/survey_optimization/simulate_survey.R): Simulates simple random, stratified random, and systematic surveys and calculates evaluation metrics. 

### 5) [Figure Production](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/figures): Create figures for use in manuscripts and presentations.

### 6) [Manuscript Production](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/manuscript): Creates manuscript for publication using RMarkdown. Don't know whether this should be shown when made public though...
