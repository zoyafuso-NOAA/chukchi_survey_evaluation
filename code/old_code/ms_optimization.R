###############################################################################
## Project:      Chukchi optimization 
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Single Species Optimization
## 
## Notes         RVersion 4.0.3
###############################################################################
rm(list = ls())

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(sp)
library(RColorBrewer)
library(raster)

##################################################
####  Load input data
##################################################
load("data/survey_opt_data/optimization_data.RData")
load("results/otter_trawl/survey_opt/SS/ss_knit_results_otter.RData")
load("results/beam_trawl_2012_2019/survey_opt/SS/ss_knit_results_beam.RData")
source("modified_functions/plot_survey_opt_map.R")
source("modified_functions/calc_expected_CV.R")
curr_dir <- getwd()

##################################################
####  
##################################################

for (igear in c("otter", "beam")) {
    
  frame_df <- get(paste0("frame_df_", igear))
  n_spp <-  get(paste0("n_spp_", igear))
  ss_strs_cvs <- get(paste0("ss_strs_cvs_", igear))
  
  ##################################################
  ####  Initialize at the simple random sample CV given 100 stations
  ##################################################
  srs_n <- 100
  srs_stats <- SamplingStrata::buildStrataDF(
    dataset = cbind( frame_df[, -grep(x = names(frame_df), pattern = "X")],
                     X1 = 1))
  ## SRS statistics
  srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:n_spp)]
  
  ##################################################
  ####  CV 
  ##################################################
  cv <- list()
  cv[["DOM"]] <- 1
  for (ispp in 1:n_spp) cv[[paste0("CV", ispp)]] <- as.numeric( srs_cv[ispp] )
  cv[["domainvalue"]] <- 1
  cv <- as.data.frame(cv)
  
  temp_n <- 0
  run <- 1
  
  ##################################################
  ####  Run optimization at a wide range of precision constraints
  ##################################################
  while (temp_n < 100) {
    
    result_dir <- paste0(curr_dir, "/results/",
                         ifelse(test = igear == "otter",
                                yes = "otter_trawl",
                                no = "beam_trawl_2012_2019"),
                         "/survey_opt/MS/Run_", run, "/")
    
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    setwd(result_dir)
    
    ## Run optimization
    par(mfrow = c(6, 6), mar = c(2, 2, 0, 0))
    solution <- SamplingStrata::optimStrata(method = "continuous",
                                            errors = cv, 
                                            framesamp = frame_df,
                                            iter = 300,
                                            pops = 100,
                                            elitism_rate = 0.1,
                                            mut_chance = 1 / (5 + 1),
                                            nStrata = 5,
                                            showPlot = T,
                                            writeFiles = T)
    
    ## Save optimization output
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
    
    plot_survey_opt_map(file_name = paste0("solution.png"),
                        grid_object =  grid_pts,
                        sol_by_cell = plot_solution, 
                        allocations = sum_stats$Allocation,
                        draw_stations = TRUE)
    
    ##################################################
    ####   Save output
    ##################################################
    result_list <- list(solution = solution,
                        sum_stats = sum_stats,
                        cvs = as.numeric(calc_expected_CV(sum_stats)),
                        n = sum(sum_stats$Allocation),
                        sol_by_cell = plot_solution)
    save(list = "result_list", file = "result_list.RData")
    
    ##################################################
    ####   Reduce CV and rerun again
    ##################################################
    cv <- list()
    cv[["DOM"]] <- 1
    for (ispp in 1:n_spp) cv[[paste0("CV", ispp)]] <- 
      result_list$cvs[ispp] * (0.9) + ss_strs_cvs[ispp] * (0.1)
    cv[["domainvalue"]] <- 1
    cv <- as.data.frame(cv)
    
    temp_n <- sum(sum_stats$Allocation)
    run <- run + 1
  }
}


