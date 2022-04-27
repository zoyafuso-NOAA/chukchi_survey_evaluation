##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Simulate Systematic Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate systematci survey designs under varying levels of 
##                      survey effort (~40 - 160 stations). 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(sp)
library(rgeos)
library(raster)
library(rgdal)
library(abind)
library(FishStatsUtils)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load input data ----
## Define areas of cells, total area
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/survey_opt_data/optimization_data.RData")
cell_area <- FishStatsUtils::chukchi_sea_grid[, "Area_in_survey_km2"]
total_area <- sum(cell_area)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create Chukchi mask ----
##   Create an outline of the Chukchi grid in aea projection
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
chukchi_mask <-
  rgdal::readOGR("data/spatial_data/survey_boundary/CHUKCHI_2012.shp")
crs(chukchi_mask) <- aea_crs
grid_pts_aea <- sp::SpatialPointsDataFrame(coords = grid_pts@coords,
                                           data = data.frame(id = 1:n_cells),
                                           proj4string = aea_crs)

##################################
## Set constants ----
## Choose the gear we are simulating. This choice will set the species list, 
## years we are simulating, and true index.
##
## For each gear, we are only simulating the most current year of data.
##################################
igear = c("otter", "beam")[1]

n_spp <-  get(paste0("n_spp_", igear))
spp_list <-  get(paste0("spp_list_", igear))
iyear <- c("beam" = "2019", "otter" = "2012")[igear]

true_index <- get(paste0("true_index_", igear))
true_index <- true_index[nrow(true_index), ]

##################################
## Simulated densities ----
## Collate simulated densities across species for a given igear.
##################################
ms_dens <- array(dim = c(n_spp, n_cells, n_iters), 
                 dimnames = list(spp_list, NULL))

for (ispp in 1:n_spp) { ## Loop over species -- start
  load(paste0("results/chukchi_", igear,
              "/vast_fits/", spp_list[ispp], "/simulated_densities/sim_data_",
              iyear, "_simtype1.RData"))
  
  ms_dens[ispp, , ] <- get(paste0("sim_data_", iyear, "_simtype1"))
} ## Loop over species -- end

##################################
## Initial result objects ----
## index_sys: estimated sample index
## cv_sys: estimated sample cv
## rb_index_sys: rel. bias of the est. index relative to the true index
## sys_settings: sample size, because we are doing a random start systematic
##               design, iterations, won't be the same exact sample size
##               for a given grid resolution. 
## Startng resultion at 65 km
temp_res <- 65000
##################################
index_sys <- cv_sys <- rb_index_sys <- list()
sys_settings <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Simulate surveys ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
set.seed(234)
while (temp_res >= 35000) {
  ## Temporary result objects for a given grid resolution
  temp_index_sys <- temp_cv_sys <- temp_rb_index_sys <-
    array(dim = c(n_spp, n_iters),
          dimnames = list(spp_list, NULL))
  for (iter in 1:n_iters) {
    
    ## Make station grid based on latlon coords, random start
    temp_grid <- sp::makegrid(x = grid_pts_aea, 
                              cellsize = temp_res, 
                              pretty = F, 
                              offset = runif(n = 2))
    temp_grid <- sp::SpatialPoints(coords = coordinates(temp_grid),
                                   proj4string = aea_crs)
    temp_grid <- raster::crop(x = temp_grid, y = chukchi_mask)
    
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
    sys_sd_LO5 <- sqrt(colSums(neighbor_var) / length(min_distances)^2)
    sys_mean <- rowMeans(ms_dens[, grid_idx, iter])
    
    ## Record survey statistics
    temp_index_sys[, iter] <- sys_mean * total_area
    temp_cv_sys[, iter] <- sys_sd_LO5 / sys_mean
    temp_rb_index_sys[, iter] <- 
      100 * (temp_index_sys[, iter] - true_index) / true_index
    
    if (iter%%50 == 0)
      print(paste0("Done with iter ", iter, 
                   ", resolution = ", temp_res, ", n = ", temp_n))
  }
  
  index_sys <- c(index_sys, list(temp_index_sys))
  cv_sys <- c(cv_sys, list(temp_cv_sys))
  rb_index_sys <- c(rb_index_sys, list(temp_rb_index_sys))
  
  print(paste("Done with resolution:", temp_res/1000, "km"))
  temp_res <- temp_res - 5000
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge results together ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
index_sys <- abind::abind(index_sys, along = 0)
cv_sys <- abind::abind(cv_sys, along = 0)
rb_index_sys <- abind::abind(rb_index_sys, along = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Performance Metrics ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
true_cv_sys <- sweep(x = apply(X = index_sys, MARGIN = 1:2, FUN = sd), 
                     MARGIN = 2, 
                     STATS = true_index, FUN = "/")

rrmse_cv_sys <- sqrt(apply(X = sweep(x = cv_sys,
                                     MARGIN = 1:2,
                                     STATS = true_cv_sys,
                                     FUN = "-") ^ 2, ## Square Error
                           MARGIN = 1:2,
                           FUN = mean) ## Mean Square Error
) /## Root Mean Square Error, Normalize by mean(cvs)
  apply(X = cv_sys, MARGIN = 1:2, FUN = mean) 

## Bias of CV relative to the True CV
bias_cv_sys <- sweep(x = cv_sys,
                     MARGIN = 1:2,
                     STATS = true_cv_sys,
                     FUN = "-")

## Relative Bias of CV relative to the True CV  
rb_cv_sys <- 100 * sweep(x = bias_cv_sys,
                         MARGIN = 1:2,
                         STATS = true_cv_sys,
                         FUN = "/")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
save(list = c("index_sys", "cv_sys", "rb_index_sys", "true_cv_sys", 
              "rrmse_cv_sys", "rb_cv_sys", "sys_settings"), 
     file = paste0("results/chukchi_", igear, "/sys_survey_sim_results.RData"))
