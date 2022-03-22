## Workflow:

### 1) [Data Inputs](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/data_processing)

* Pull Chukchi 1990 and 2012 otter (83-112 eastern bottom trawl) catch and effort data from RACE database.
* Pull Chukchi 2012, 2017, and 2019 beam trawl data from IERL survey.

### 2) [Operating Model](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/VAST_fitting)

* Fit single-species VAST models with different model configurations (e.g., spatial/spatiotemporal fields, observation models). Save the best (lowest-aic) model, run diagnostics, and simulate densities.  

### 3) [Stratified Survey Optimization](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/survey_optimization)

* Conduct stratified-random design optimization.

### 4) Survey Design Evaluation

* Simulate different survey designs and calculate performance metrics.

### 5) [Figure Production](https://github.com/zoyafuso-NOAA/Arctic_GF_OM/tree/main/code/figures)

* Create figures for use in manuscripts and presentations.

### 6) Manuscript Production

* Create manuscript for publication. 
