## Workflow:

### 1) [Data Inputs](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/tree/main/code/00_data_processing)

* [Process_AK_BTS_data_arctic.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/00_data_processing/Process_AK_BTS_data_arctic.R): Pulls Chukchi 1990 and 2012 otter (83-112 eastern bottom trawl) catch and effort data from RACE database.
* [Process_AK_IERL_data_arctic.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/00_data_processing/Process_AK_IERL_data_arctic.R): Pulls Chukchi 2012, 2017, and 2019 beam trawl data from IERL survey.
* [aggregates_species.csv](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/data/fish_data/aggregates_species.csv): table defining species taxon groups, called in the two above scripts. 

### 2) [Analysis](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/tree/main/code/01_analysis)

* [VAST_fitting_general.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/01_analysis/VAST_fitting_general.R): Fits single-species VAST models with different model configurations (e.g., spatial/spatiotemporal fields, observation models). Saves the best (lowest-aic) model, runs diagnostics, and simulates densities.  
* [optimization_data.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/01_analysis/optimization_data.R): Sets up constants used in the (subsequent) STRS design optimization and survey simulation scripts.
* [combined_ss_ms_optimization.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/01_analysis/combined_ss_ms_optimization.R): Conducts stratified-random design optimization for both gears.
* [simulate_survey.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/code/01_analysis/simulate_survey.R): Simulates simple random, stratified random, and systematic surveys and calculates evaluation metrics. 

### 3) [Figure Production](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/tree/main/code/02_figures): Create figures for use in manuscripts and presentations.

### 4) [Manuscript Production](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/tree/main/manuscript): Creates manuscript for publication using RMarkdown. 

* [run.R](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/manuscript/run.R): main script that produces the manuscript documents, creates appendices. The main use of this script is in creating and saving different versions of the manuscript (e.g., Version_0, Version_1, etc.) as the manuscript updates throughout the review process. 
* [Oyafuso_etal_Chukchi_MS.Rmd](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/manuscript/Oyafuso_etal_Chukchi_MS.Rmd): RMarkdown script for the manuscript content. 
* [Oyafuso_etal_Chukchi_survey_testing.bib](https://github.com/zoyafuso-NOAA/chukchi_survey_evaluation/blob/main/manuscript/Oyafuso_etal_Chukchi_survey_testing.bib): .bib file containing the references used in the manuscript. 
