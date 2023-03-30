##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Simulate Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate different survey designs under varying levels of 
##                      survey effort.
##
##                Survey Designs: 1) Simple Random Desgin
##                                2) Fixed Grid Systematic Design
##                                3) Stratified Random Design
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
library(terra)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load input data ----
## Define areas of cells, total area
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/survey_opt_data/optimization_data_nbs.RData")
cell_area <- FishStatsUtils::northern_bering_sea_grid[, "Area_in_survey_km2"]
total_area <- sum(cell_area)

grid_pts <- terra::vect("data/survey_opt_data/grid_pts_nbs.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create an outline of the Chukchi grid in aea projection. This mask is used
##       when creating the systematic grids. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbs_mask <-
  rgdal::readOGR("data/spatial_data/survey_boundary/NBS_2019.shp")
crs(nbs_mask) <- aea_crs

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set constants ----
## Choose the gear we are simulating. This choice will set the species list, 
## years we are simulating, and true index.
##
## For each gear, we are only simulating the most current year of data.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

igear = c("otter", "beam")[1]

true_index <- true_index[nrow(true_index), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulated densities ----
## Collate simulated densities across species for a given igear.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
ms_dens <- array(dim = c(n_spp, n_cells, n_iters), 
                 dimnames = list(spp_list, NULL))

for (ispp in 1:n_spp) { ## Loop over species -- start
  load(paste0("results/nbs_", igear,
              "/vast_fits/", spp_list$taxon[ispp], "/simulated_densities/sim_data_",
              iyear, "_simtype1.RData"))
  
  ms_dens[ispp, , ] <- get(paste0("sim_data_", iyear, "_simtype1"))
} ## Loop over species -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate SRS ----
## Simulate Simple Random Designs at varying sampling efforts   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
target_n <- seq(from = 55, to = 175, by = 15)
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

rm(srs_tau, srs_sd, srs_sample, srs_sample_idx, isample, iter)

##################################
##   Simulate SYS surveys ----
## index_sys: estimated sample index
## cv_sys: estimated sample cv
## rb_index_sys: rel. bias of the est. index relative to the true index
## sys_settings: sample size, because we are doing a random start systematic
##               design, iterations, won't be the same exact sample size
##               for a given grid resolution. 
## Startng resultion at 65 km
temp_res <- 65000
##################################
index_sys <- cv_sys <- cv_sys_LO5 <- rb_index_sys <- list()
sys_settings <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create points within NBS mask
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
grid_pts_aea <- sp::SpatialPointsDataFrame(coords = grid_pts@coords,
                                           data = data.frame(id = 1:n_cells),
                                           proj4string = aea_crs)


set.seed(54904)
while (temp_res >= 35000) {
  
  ## Temporary result objects for a given grid resolution
  temp_index_sys <- temp_cv_sys <- temp_cv_sys_LO5 <- temp_rb_index_sys <-
    array(dim = c(n_spp, n_iters),
          dimnames = list(spp_list, NULL))
  
  ## Make station grid based on latlon coords, random start
  temp_grid <- sp::makegrid(x = grid_pts_aea, 
                            cellsize = temp_res, 
                            pretty = F)
  temp_grid <- sp::SpatialPoints(coords = coordinates(temp_grid),
                                 proj4string = aea_crs)
  temp_grid <- terra::intersect(x = temp_grid, y = chukchi_mask)
  
  plot(temp_grid)
  
  ## Calculate sample size
  temp_n <- length(temp_grid)
  sys_settings <- rbind(sys_settings,
                        data.frame(n = temp_n,
                                   res = temp_res))
  
  ## Calculate the grid cell that falls under the each station
  grid_idx <- apply(X = rgeos::gDistance(spgeom1 = grid_pts_aea,
                                         spgeom2 = temp_grid,
                                         byid = TRUE),
                    MARGIN = 1,
                    FUN = which.min)
  
  ## Calculate the nearest neighbor stations for each station
  distances <- as.matrix(dist(grid_pts_aea[grid_idx, ]@coords))
  min_distances <- apply(X = distances,
                         MARGIN = 1,
                         FUN = function(x) {
                           y = round(x[x > 0])
                           return(names(which( y == min(y) )))
                         })
  
  neighbor_var <- matrix(nrow = length(min_distances), ncol = n_spp,
                         dimnames = list(NULL, spp_list))
  
  for (iter in 1:n_iters) {
    
    for (icell in 1:length(min_distances) ){  ## Loop over stations -- start
      
      ## Pull the station indices of the station and its neighbors
      idx <- as.integer(c(icell, min_distances[[icell]]))
      
      ## Pull the location of the idx on the grid from sys_idx and get the
      ## simulated densities
      neighbors <- ms_dens[, grid_idx[idx], iter]
      
      ## Calculate the neighborhood variance of the station
      neighbor_var[icell, ] <-  apply(X = neighbors,
                                      MARGIN = 1,
                                      FUN = function(x) mean((x - mean(x))^2))
    } ## Loop over stations -- end
    
    ## Calculate survey statistics
    sys_sample <- ms_dens[, grid_idx, iter]
    sys_mean <- rowMeans(ms_dens[, grid_idx, iter])
    sys_sd_LO5 <- sqrt(colSums(neighbor_var) / length(min_distances)^2)
    sys_index <- rowMeans(sys_sample) * total_area
    sys_sd <- sqrt(apply(X = sys_sample, MARGIN = 1, FUN = var) *
                     total_area^2 / temp_n)
    
    ## Record survey statistics
    temp_index_sys[, iter] <- sys_index
    temp_cv_sys[, iter] <- sys_sd / sys_index
    temp_cv_sys_LO5[, iter] <- sys_sd_LO5 / sys_mean
    temp_rb_index_sys[, iter] <- 
      100 * (temp_index_sys[, iter] - true_index) / true_index
    
    if (iter%%50 == 0)
      print(paste0("Done with iter ", iter, 
                   ", resolution = ", temp_res, ", n = ", temp_n))
  }
  
  index_sys <- c(index_sys, list(temp_index_sys))
  cv_sys <- c(cv_sys, list(temp_cv_sys))
  cv_sys_LO5 <- c(cv_sys_LO5, list(temp_cv_sys_LO5))
  rb_index_sys <- c(rb_index_sys, list(temp_rb_index_sys))
  
  print(paste("Done with resolution:", temp_res/1000, "km"))
  temp_res <- temp_res - 2500
}

index_sys <- abind::abind(index_sys, along = 0)
rb_index_sys <- abind::abind(rb_index_sys, along = 0)
cv_sys <- abind::abind(cv_sys, along = 0)
cv_sys_LO5 <- abind::abind(cv_sys_LO5, along = 0)

rm(sys_sample, sys_mean, sys_sd_LO5, sys_index, sys_sd, 
   temp_grid, distances, 
   temp_res, temp_rb_index_sys, temp_n,
   temp_index_sys, temp_cv_sys, temp_cv_LO5, temp_cv_sys_LO5)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Simumulate STRS ----
##   5-strata multispecies STRS optimization
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
index_ms_strs <- cv_ms_strs <- rb_index_ms_strs <-
  array(dim = c(length(target_n), n_spp, n_iters),
        dimnames = list(paste0("ms_strs_n = ", target_n), spp_list, NULL))
n_strata <- c("otter" = 3, "beam" = 4)[igear]

load(paste0("results/chukchi_", igear,
            "/survey_opt/Str_", n_strata, "/result_list.RData"))
load(paste0("results/chukchi_", igear,
            "/survey_opt/Str_", n_strata, "/allocations.RData"))

set.seed(5231)
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
##   Performance Metrics ----
##   True CV: sd of the 1000 estimated index simulation replicates divided by
##            the true index
##   RRMSE of CV: 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (idesign in c("srs", "ms_strs", "sys")) { ## Loop over design type -- start
  
  ## Pull the simulated indices from a particular idesign
  temp_index <- get(paste0("index_", idesign))
  
  ## Remove positive outliers > 3 SDs from the mean for a given combination
  ##     of species/sample size.
  temp_index <- apply(X = temp_index, 
                      MARGIN = 1:2, 
                      FUN = function(x) ifelse(test = x > mean(x) + 3*sd(x), 
                                               yes = NA, 
                                               no = x ))
  
  ## Pull the simulated cvs from a particular idesign
  temp_cv <- get(paste0("cv_", idesign))
  if (idesign == "sys") temp_cv_LO5 <- cv_sys_LO5
  
  ## Remove positive outliers > 3 SDs from the mean for a given combination
  ##     of species/sample size.
  temp_cv <- apply(X = temp_cv, 
                   MARGIN = 1:2, 
                   FUN = function(x) ifelse(test = x > mean(x) + 3*sd(x), 
                                            yes = NA, 
                                            no = x ))
  
  if (idesign == "sys") 
    temp_cv_LO5 <- apply(X = temp_cv_LO5, 
                         MARGIN = 1:2, 
                         FUN = function(x) ifelse(test = x > mean(x) + 3*sd(x), 
                                                  yes = NA, 
                                                  no = x ))
  ## True CV
  temp_true_cv <- sweep(x = apply(X = temp_index, 
                                  MARGIN = 2:3, 
                                  FUN = sd, 
                                  na.rm = TRUE),
                        MARGIN = 2,
                        STATS = true_index,
                        FUN = "/")
  
  ## RRMSE of CV
  temp_rrmse_cv <- sqrt(apply(X = sweep(x = temp_cv,
                                        MARGIN = 2:3,
                                        STATS = temp_true_cv,
                                        FUN = "-") ^ 2, ## Square Error
                              MARGIN = 2:3,
                              FUN = mean,
                              na.rm = TRUE) ## Mean Square Error
  ) /## Root Mean Square Error
    apply(X = temp_cv, 
          MARGIN = 2:3, 
          FUN = mean, 
          na.rm = TRUE) ## Normalize by mean(cvs)
  
  if (idesign == "sys") 
    temp_rrmse_cv_LO5 <- sqrt(apply(X = sweep(x = temp_cv_LO5,
                                              MARGIN = 2:3,
                                              STATS = temp_true_cv,
                                              FUN = "-") ^ 2, ## Square Error
                                    MARGIN = 2:3,
                                    FUN = mean,
                                    na.rm = TRUE) ## Mean Square Error
    ) /## Root Mean Square Error
    apply(X = temp_cv_LO5, 
          MARGIN = 2:3, 
          FUN = mean, 
          na.rm = TRUE) ## Normalize by mean(cvs)
  
  ## Bias of CV relative to the True CV
  temp_bias_cv <- sweep(x = temp_cv,
                        MARGIN = 2:3,
                        STATS = temp_true_cv,
                        FUN = "-")
  
  ## Relative Bias of CV relative to the True CV  
  temp_rb_cv <- 100 * sweep(x = temp_bias_cv,
                            MARGIN = 2:3,
                            STATS = temp_true_cv,
                            FUN = "/")
  
  if (idesign == "sys") {
    
    ## Bias of CV relative to the True CV
    temp_bias_cv_LO5 <- sweep(x = temp_cv_LO5,
                              MARGIN = 2:3,
                              STATS = temp_true_cv,
                              FUN = "-")
    
    ## Relative Bias of CV relative to the True CV  
    temp_rb_cv_LO5 <- 100 * sweep(x = temp_bias_cv_LO5,
                                  MARGIN = 2:3,
                                  STATS = temp_true_cv,
                                  FUN = "/")
  }
  
  ## Assign temporary variables to unique names consisting of the metric with 
  ## the design type.
  for (imetric in c("true_cv", 
                    "rrmse_cv", 
                    "rb_cv")) { ## Loop over metric -- start
    assign(x = paste0(imetric, "_", idesign), 
           value = get(paste0("temp_", imetric)) )
    
    rm(list = paste0("temp_", imetric))
    
  } ## Loop over metric -- end
  
  if (idesign == "sys") {
    rrmse_cv_sys_LO5 <- temp_rrmse_cv_LO5
    rb_cv_sys_LO5 <- temp_rb_cv_LO5
    
    rm(temp_rrmse_cv_LO5, temp_bias_cv_LO5, temp_rb_cv_LO5)
  }
} ## Loop over design type -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(file = paste0("results/chukchi_", igear, 
                   "/random_survey_sim_results.RData"),
     list = c(as.vector(sapply(X = c("true_cv_", "rrmse_cv_", 
                                     "index_", "cv_", 
                                     "rb_index_", "rb_cv_"), 
                               FUN = function(x) 
                                 paste0(x, c("srs", "ms_strs", "sys")))),
              "target_n", "sys_settings", 
              "cv_sys_LO5", "rb_cv_sys_LO5", "rrmse_cv_sys_LO5"))
