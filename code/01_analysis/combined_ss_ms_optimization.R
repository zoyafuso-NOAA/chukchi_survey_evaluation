##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:      Chukchi optimization 
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Optimization scheme for a multispecies stratified design
## 
## Notes         RVersion 4.0.2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages ----
##  Install a forked version of the SamplingStrata Package from 
##  zoyafuso-NOAA's Github page
##
##  Import other required packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(sp)
library(RColorBrewer)
library(raster)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load Input Data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
load("data/survey_opt_data/optimization_data.RData")
source("modified_functions/plot_survey_opt_map.R")
source("modified_functions/calc_expected_CV.R")
curr_dir <- getwd()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Assign Gear ----
##  Define species list and years given the gear
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
iregion = "chukchi"
igear = "otter"
n_spp <- get(paste0("n_spp_", igear))
spp_list <- get(paste0("spp_list_", igear))
frame_df <- get(paste0("frame_df_", igear))
n_years <- c("otter" = 2, "beam" = 3)[igear]

stratas <- 3:4
total_n <- seq(from = 55, to = 200, by = 15)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Conduct Optimization ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

for (istrata in stratas[1]) {
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   1) Initialize Optimization ----
  ##   Calculate spatiotemporal variance given a simple random sample with 
  ##   200 stations for each species. These variances are used as the starting
  ##   points for the survey design optimization because we assume the optimized
  ##   stratified design will produce better variances than the simple random
  ##   design for a given level of sampling effort.
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  srs_n <- 200
  srs_stats <- SamplingStrata::buildStrataDF(
    dataset = cbind( frame_df[, -grep(x = names(frame_df), pattern = "X")],
                     X1 = 1))
  ## SRS statistics
  srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:n_spp)]
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   2) Create CV Input ----
  ##   CV dataframe for optimization input 
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  cv <- list()
  cv[["DOM"]] <- 1
  for (ispp in 1:n_spp) cv[[paste0("CV", ispp)]] <- as.numeric( srs_cv[ispp] )
  cv[["domainvalue"]] <- 1
  cv <- as.data.frame(cv)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   3) Create result directory ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  result_dir <- paste0(curr_dir, "/results/",
                       iregion, "_", igear, "/survey_opt",
                       "/Str_", istrata, "/")
  
  if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
  setwd(result_dir)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   4) Run initial optimization ----
  ##   subject to the SRS CVs
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  par(mfrow = c(6, 6), mar = c(2, 2, 0, 0))
  solution <- SamplingStrata::optimStrata(method = "continuous",
                                          errors = cv, 
                                          framesamp = frame_df,
                                          iter = 200,
                                          pops = 100,
                                          elitism_rate = 0.1,
                                          mut_chance = 1 / (istrata + 1),
                                          nStrata = istrata,
                                          showPlot = T,
                                          writeFiles = T)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   5) Clean up results ---- 
  ##   Clean optimization output and save
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
  solution$aggr_strata <- 
    solution$aggr_strata[order(solution$aggr_strata$DOM1,
                               solution$aggr_strata$STRATO), ]
  
  sum_stats <- SamplingStrata::summaryStrata(solution$framenew,
                                             solution$aggr_strata,
                                             progress=FALSE)
  sum_stats$stratum_id <- 1:nrow(sum_stats)
  sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
  sum_stats$Wh <- sum_stats$Population / nrow(frame_df)
  sum_stats <- cbind(sum_stats,
                     subset(x = solution$aggr_strata,
                            select = -c(STRATO, N, COST, CENS, DOM1, X1)))
  
  plot_solution <- solution$indices$X1
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   6) Save results ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  result_list <- list(solution = solution,
                      sum_stats = sum_stats,
                      cvs = as.numeric(calc_expected_CV(sum_stats)),
                      n = sum(sum_stats$Allocation),
                      sol_by_cell = plot_solution)
  save(list = "result_list", file = "result_list.RData")
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   7) Single-Species Optimization ----
  ##   Calculate single-species CV subject to the initial stratification
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  ss_sample_allocations <- expand.grid(n = total_n, species = spp_list)
  for (ispp in 1:n_spp) {
    temp_n <- result_list$n
    
    ## Subset density data for species ispp
    ss_df <- subset(x = frame_df, 
                    select = c("domainvalue", "id", "WEIGHT", "X1", "X2", 
                               paste0("Y", ispp), paste0("Y", ispp, "_SQ_SUM")))
    names(ss_df)[grep(x = names(ss_df), pattern = "Y")] <- c("Y1", "Y1_SQ_SUM")
    
    for (isample in total_n){
      
      ## Create CV inputs to the Bethel algorithm; initialize at SRS CV
      error_df <- data.frame("DOM" = "DOM1",
                             as.numeric(srs_cv)[ispp],
                             "domainvalue"  = 1)
      names(error_df)[2] <- "CV1"
      
      ## subset stratum stats for the species of interest as inputs to the 
      ## Bethel algorithm
      temp_stratif <- 
        solution$aggr_strata[, c("STRATO", "N", 
                                 paste0("M", ispp), paste0("S", ispp), 
                                 "COST", "CENS", "DOM1", "X1" , "SOLUZ"
        )]
      temp_stratif$N <- temp_stratif$N / n_years
      temp_stratif$DOM1 <- 1
      names(temp_stratif)[3:4] <- paste0(c("M", "S"), 1)
      
      ## run bethel at the SRS CV
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif, 
        realAllocation = T, 
        printa = T)
      
      ## Save the current n and cv constraint
      temp_n <- sum(ceiling(temp_bethel))
      updated_cv_constraint <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
      
      ## modify CV, rerun bethel, save temp_n, run until temp_n == isample
      while (temp_n != isample){
        over_under <- temp_n > isample
        CV_adj <- ifelse(over_under == TRUE, 
                         yes = 1.01,
                         no = 0.999)
        
        updated_cv_constraint <- updated_cv_constraint * CV_adj
        
        error_df[, "CV1"] <- updated_cv_constraint
        
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df, 
                                              printa = TRUE)
        
        temp_n <- sum(as.numeric(temp_bethel))
        updated_cv_constraint <- 
          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
        print(paste0("n = ", temp_n, ", ", updated_cv_constraint) )
      }
      
      ## Save the CV and station allocations that corresponds to isample
      temp_idx <- ss_sample_allocations$n == isample & 
        ss_sample_allocations$species == spp_list[ispp]
      
      ss_sample_allocations[temp_idx, "CV"] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
      
      ss_sample_allocations[temp_idx, paste0("Str_", 1:length(temp_bethel))] <- 
        as.integer(temp_bethel)
      
    }
  }
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   8) Adjust MS solution ----
  ##   Optimize allocation across a range of sample sizes, given the original
  ##   stratification.
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  ms_sample_allocations <- expand.grid(n = total_n)
  temp_n <- result_list$n
  
  for (isample in total_n){
    
    ## Subset lower limits of CVs from the ss cvs
    ss_cvs <- subset(ss_sample_allocations, n == isample)$CV
    
    ## CV dataframe input to Bethel algorithm
    error_df <-  data.frame("DOM" = "DOM1",
                            srs_cv,
                            "domainvalue"  = 1)
    names(error_df)[2:(1 + n_spp)] <- paste0("CV", 1:n_spp)
    
    ## Stratum statistics input to the Bethel algorithm
    temp_stratif <- solution$aggr_strata
    temp_stratif$N <- temp_stratif$N / n_years
    temp_stratif$DOM1 <- 1
    
    ## Run Bethel algorithm and save current n and cv constraints
    temp_bethel <- SamplingStrata::bethel(
      errors = error_df,
      stratif = temp_stratif, 
      realAllocation = T, 
      printa = T)
    
    temp_n <- sum(ceiling(temp_bethel))
    updated_cv_constraint <- 
      as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
    
    ## Rerun Bethel algorithm, modifying the CVs relative to the distances
    ## between the SRS and SS CVs given isample stations. First we calculate 
    ## CVs calculated under SRS for each species with isample stations
    temp_srs_var <- 
      srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - isample / n_cells) / isample
    temp_srs_cv <- sqrt(temp_srs_var) / srs_stats[, paste0("M", 1:n_spp)]
    
    while (temp_n != isample) {
      over_under <- temp_n > isample
      
      ## If the current n is < isample, decrease the CV by a small amount
      ## relative to the distance between the current CV and the SS CV
      if (over_under == FALSE) {
        CV_adj <- 0.95
        updated_cv_constraint <- 
          updated_cv_constraint * (CV_adj) + ss_cvs * (1  - CV_adj)
      }
      
      ## If the current n is > isample, increase the CV by a small amount
      ## relative to the distance between the current CV and the SRS CV
      if(over_under == TRUE) {
        CV_adj = .05
        updated_cv_constraint <- 
          temp_srs_cv * (CV_adj) + updated_cv_constraint * (1  - CV_adj)
      }
      
      ## Update the CV dataframe input with the updated_cv_constraint
      error_df[, paste0("CV", 1:n_spp)] <- updated_cv_constraint
      
      ## Rerun Bethel algorithm
      temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                            errors = error_df, 
                                            printa = TRUE)
      
      ## Save sample size and CV constraint
      temp_n <- sum(as.numeric(temp_bethel))
      updated_cv_constraint <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
      
      ## Print out result to console
      print(paste0("n = ", temp_n) )
    }
    
    ## Save optimized CV 
    temp_idx <- ms_sample_allocations$n == isample
    
    ms_sample_allocations[temp_idx, paste0("CV", 1:n_spp)] <- 
      as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
    
    ms_sample_allocations[temp_idx, paste0("Str_", 1:length(temp_bethel))] <- 
      as.integer(temp_bethel)
    
  }
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   9) Save solution image ----
  ##   Save an image of the solution with 200 stations
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  plot_survey_opt_map(file_name = paste0("solution.png"),
                      grid_object =  grid_pts,
                      sol_by_cell = plot_solution, 
                      allocations = as.numeric(ms_sample_allocations[ms_sample_allocations$n == 190, 
                                                                     paste0("Str_", 1:length(temp_bethel))]),
                      draw_stations = TRUE)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   10) Save allocations ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  save(list = c("ss_sample_allocations", "ms_sample_allocations"), 
       file = "allocations.RData")
}

