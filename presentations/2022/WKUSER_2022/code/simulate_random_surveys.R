##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Simulate Random Desgin Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate different survey designs under varying levels of 
##                      survey effort.
##
##                Survey Designs: Simple Random, Stratified Random
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(sp)
library(rgeos)
library(raster)
library(rgdal)
library(viridis )
library(FishStatsUtils)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load input data ----
## Define areas of cells, total area
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/survey_opt_data/optimization_data.RData")
cell_area <- FishStatsUtils::chukchi_sea_grid[, "Area_in_survey_km2"]
total_area <- sum(cell_area)


load(paste0("presentations/2022/WKUSER_2022/results/Str_3/result_list.RData"))
load(paste0("presentations/2022/WKUSER_2022/results/Str_3/allocations.RData"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create an outline of the Chukchi grid in aea projection
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chukchi_mask <-
  rgdal::readOGR("data/spatial_data/survey_boundary/CHUKCHI_2012.shp")
crs(chukchi_mask) <- aea_crs

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set constants ----
## Choose the gear we are simulating. This choice will set the species list, 
## years we are simulating, and true index.
##
## For each gear, we are only simulating the most current year of data.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
igear = "otter"
spp_list <- get(paste0("spp_list_", igear))
spp_list <- spp_list[spp_list %in% c("Arctic cod", "snow crab",
                                     "Bering flounder", "saffron cod", 
                                     "Alaska plaice", "yellowfin sole")]
n_spp <- length(spp_list)
iyear <- c("beam" = "2019", "otter" = "2012")[igear]

true_index <- get(paste0("true_index_", igear))
true_index <- true_index[nrow(true_index), names(spp_list)]
names(true_index) <- spp_list

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Predicted densities ----
## Collate predicted densities across species for a given igear.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
ms_dens <- array(dim = c(n_spp, n_cells, n_iters), 
                 dimnames = list(spp_list, NULL, NULL))

for (ispp in 1:n_spp) { ## Loop over species -- start
  # load(paste0("results/chukchi_", igear,
  #             "/vast_fits/", spp_list[ispp], "/fit.RData"))
  # ms_dens[ispp, ] <- fit$Report$D_gct[, 1, "2012"]
  load(paste0("results/chukchi_", igear,
              "/vast_fits/", spp_list[ispp], "/simulated_densities/sim_data_",
              iyear, "_simtype1.RData"))
  ms_dens[ispp, , ] <- get(paste0("sim_data_", iyear, "_simtype1"))

} ## Loop over species -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate SRS ----
## Simulate Simple Random Designs at varying sampling efforts   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
target_n <- seq(from = 50, to = 200, by = 50)
index_srs <- cv_srs <- rb_index_srs <-
  array(dim = c(length(target_n), n_spp, n_iters),
        dimnames = list(paste0("srs_n = ", target_n), spp_list, NULL))

set.seed(234)
for (isample in 1:length(target_n)) { ## Loop over total effort -- start
  for (iter in 1:n_iters) { ## Loop over iterations -- start
    
    ## Randomly sample grid indices
    srs_sample_idx <- sample(x = 1:n_cells,
                             size = target_n[isample])
    srs_sample <- ms_dens[, srs_sample_idx, iter]
    
    ## Calculate survey statistics
    srs_tau <- total_area * rowMeans(srs_sample)
    srs_sd <- sqrt(apply(X = srs_sample, MARGIN = 1, FUN = var) *
                     total_area^2 / target_n[isample])
    
    ## Record Survey Performance Metrics
    index_srs[isample, , iter] <- srs_tau 
    rb_index_srs[isample, , iter] <- 
      100 * (srs_tau - true_index) / true_index
    cv_srs[isample, , iter] <-  srs_sd / srs_tau
    
  }  ## Loop over iterations -- end
} ## Loop over total effort -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Simumulate STRS ----
##   3-strata multispecies STRS optimization
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
index_ms_strs <- cv_ms_strs <- rb_index_ms_strs <-
  array(dim = c(length(target_n), n_spp, n_iters),
        dimnames = list(paste0("ms_strs_n = ", target_n), spp_list, NULL))
n_strata <- 3

set.seed(231)
for (itarget in target_n){ ## Loop over sampling effort -- start
  
  ## STRS settings
  solution <- result_list$sol_by_cell
  allocation <- as.integer(subset(ms_sample_allocations, 
                                  n == itarget)[, paste0("Str_", 1:n_strata)])
  Ah <- tapply(X = cell_area, INDEX = solution, FUN = sum)
  n_strs <- sum(allocation)
  
  for (iter in 1:n_iters) { ## Loop over iterations -- start
    sample_vec <- stratum_order <- c()
    for (istrata in 1:length(allocation)) { ## Loop over strata -- start
      sample_vec <- c(sample_vec,
                      sample(x = which(solution == istrata),
                             size = allocation[istrata]) )
      stratum_order <- c(stratum_order,
                         rep(x = istrata, times = allocation[istrata]))
    } ## Loop over strata -- end
    
    # Subset sub_df by which cells were chosen
    sample_df <- ms_dens[, sample_vec, iter]
    
    # Calculate STRS mean density
    strata_mean <- apply(X = sample_df, 
                         MARGIN = 1, 
                         FUN = function(x) tapply(X = x, 
                                                  INDEX = stratum_order,
                                                  FUN = mean))
    
    tau_strs <- as.numeric(t(strata_mean) %*% matrix(Ah) )
    index_ms_strs[paste0("ms_strs_n = ", itarget), , iter] <- tau_strs
    
    # Calculate STRS variance of mean density
    strata_var <- apply(X = sample_df, 
                        MARGIN = 1, 
                        FUN = function(x) 
                          tapply(X = x, 
                                 INDEX = stratum_order,
                                 FUN = function(xx)
                                   (var(xx) / (length(xx)) )) )
    
    STRS_var <-  as.numeric(t(strata_var) %*% matrix(Ah^2))
    
    # Save mean and cv of estimates across species
    cv_ms_strs[paste0("ms_strs_n = ", itarget), , iter]  <- 
      sqrt(STRS_var) / tau_strs
    
    rb_index_ms_strs[paste0("ms_strs_n = ", itarget), , iter] <- 
      100 * (tau_strs - true_index) / true_index
    
  } ## Loop over iterations -- end
  print(paste("Done with n =", itarget))
} ## Loop over sampling effort -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Simumulate STRS ----
##   3-strata single species STRS optimization
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
index_ss_strs <- cv_ss_strs <- rb_index_ss_strs <-
  array(dim = c(length(target_n), n_spp, n_iters),
        dimnames = list(paste0("ss_strs_n = ", target_n), spp_list, NULL))
n_strata <- 3

set.seed(231)
for (itarget in target_n){ ## Loop over sampling effort -- start
  
  ## STRS settings
  solution <- result_list$sol_by_cell
  Ah <- tapply(X = cell_area, INDEX = solution, FUN = sum)
  
  for (ispp in spp_list) {
    allocation <- as.integer(subset(ss_sample_allocations, 
                                    n == itarget & species == ispp)[, paste0("Str_", 1:n_strata)])
    n_strs <- sum(allocation)
    
    for (iter in 1:n_iters) { ## Loop over iterations -- start
      sample_vec <- stratum_order <- c()
      for (istrata in 1:length(allocation)) { ## Loop over strata -- start
        sample_vec <- c(sample_vec,
                        sample(x = which(solution == istrata),
                               size = allocation[istrata]) )
        stratum_order <- c(stratum_order,
                           rep(x = istrata, times = allocation[istrata]))
      } ## Loop over strata -- end
      
      # Subset sub_df by which cells were chosen
      sample_df <- ms_dens[ispp, sample_vec, iter]
      
      # Calculate STRS mean density
      strata_mean <- tapply(X = sample_df, 
                            INDEX = stratum_order,
                            FUN = mean)
      
      tau_strs <- as.numeric(t(strata_mean) %*% matrix(Ah) )
      index_ss_strs[paste0("ss_strs_n = ", itarget), ispp, iter] <- tau_strs
      
      # Calculate STRS variance of mean density
      strata_var <- tapply(X = sample_df, 
                           INDEX = stratum_order,
                           FUN = function(xx)
                             (var(xx) / (length(xx)) )) 
      
      STRS_var <-  as.numeric(t(strata_var) %*% matrix(Ah^2))
      
      # Save mean and cv of estimates across species
      cv_ss_strs[paste0("ss_strs_n = ", itarget), ispp, iter]  <- 
        sqrt(STRS_var) / tau_strs
      
      rb_index_ss_strs[paste0("ss_strs_n = ", itarget), ispp, iter] <- 
        100 * (tau_strs - true_index[ispp]) / true_index[ispp]
      
    } ## Loop over iterations -- end
    print(paste("Done with n =", itarget, ispp))
  }
  
  
  
} ## Loop over sampling effort -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Performance Metrics ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (idesign in c("srs", "ms_strs", "ss_strs")) { ## Loop over design type -- start
  
  temp_index <- get(paste0("index_", idesign))
  temp_cv <- get(paste0("cv_", idesign))
  
  ## True CV
  temp_true_cv <- sweep(x = apply(X = temp_index, MARGIN = 1:2, FUN = sd),
                        MARGIN = 2,
                        STATS = true_index,
                        FUN = "/")
  
  ## RRMSE of CV
  temp_rrmse_cv <- sqrt(apply(X = sweep(x = temp_cv,
                                        MARGIN = 1:2,
                                        STATS = temp_true_cv,
                                        FUN = "-") ^ 2, ## Square Error
                              MARGIN = 1:2,
                              FUN = mean) ## Mean Square Error
  ) /## Root Mean Square Error
    apply(X = temp_cv, MARGIN = 1:2, FUN = mean) ## Normalize by mean(cvs)
  
  ## CV of CV
  temp_cv_cv <- apply(X = temp_cv,
                      MARGIN = 1:2,
                      FUN = function(x) sd(x) / mean(x) )
  
  ## Bias of CV relative to the True CV
  temp_bias_cv <- sweep(x = temp_cv,
                        MARGIN = 1:2,
                        STATS = temp_true_cv,
                        FUN = "-")
  
  ## Relative Bias of CV relative to the True CV  
  temp_rb_cv <- 100 * sweep(x = temp_bias_cv,
                            MARGIN = 1:2,
                            STATS = temp_true_cv,
                            FUN = "/")
  
  ## Assign temporary variables to unique names consisting of the metric with 
  ## the design type.
  for (imetric in c("true_cv", "rrmse_cv", 
                    "cv_cv", "rb_cv")) { ## Loop over metric -- start
    assign(x = paste0(imetric, "_", idesign), 
           value = get(paste0("temp_", imetric)) )
    rm(list = paste0("temp_", imetric))
  } ## Loop over metric -- end
  
} ## Loop over design type -- end

rrmse_cv_srs
rrmse_cv_ms_strs


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(file = paste0("presentations/2022/WKUSER_2022/results/Str_3/",
                   "random_survey_sim_results.RData"),
     list = c(as.vector(sapply(X = c("true_cv_", "rrmse_cv_",
                                     "index_", "cv_",
                                     "rb_index_", "rb_cv_"),
                               FUN = function(x)
                                 paste0(x, c("srs", "ss_strs", "ms_strs")))),
              "target_n"))
