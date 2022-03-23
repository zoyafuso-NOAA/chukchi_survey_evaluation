###############################################################################
## Project:       Simulate Surveys Function
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate different survey designs under varying levels of 
##                      survey effort (40 - 100 stations). 
##
##                Survey Designs: Simple Random, Systematic, Stratified Random
###############################################################################
rm(list = ls())

##################################
## Import Libraries
##################################
library(sp)
library(rgeos)
library(raster)
library(rgdal)
library(viridis )
library(VAST)

##################################
## Load input data and single-species survey outputs
## Load Chukchi grid
##################################
load("data/survey_opt_data/optimization_data.RData")
# chukchi_grid <- read.csv(paste0("data/spatial_data/",
#                                 "BS_Chukchi_extrapolation_grids/",
#                                 "ChukchiThorsonGrid.csv"))

##################################################
####   convert shape area (m2) to km2, assign some area constants
##################################################  
# chukchi_grid$Area_km2 <- chukchi_grid$Shape_Area / 1000 / 1000
cell_area <- chukchi_sea_grid[, "Area_in_survey_km2"]
total_area <- sum(cell_area)

##################################################
####   Create an outline of the Chukchi grid in aea projection
##################################################  
chukchi_mask <-
  rgdal::readOGR("data/spatial_data/survey_boundary/CHUKCHI_2012.shp")
crs(chukchi_mask) <- aea_crs

##################################
## Calculate locations of systematic designs with different sampling efforts
##################################
{
  grid_pts_aea <- sp::SpatialPointsDataFrame(coords = grid_pts@coords,
                                             data = data.frame(id = 1:n_cells),
                                             proj4string = aea_crs)
  temp_n <- 0
  temp_res <- 65000
  sys_settings <- NULL
  sys_idx <- list()
  
  # png(filename = "presentations/results_11_30_2021/sys_survey_grids.png",
  #     width = 4, height = 6, units = "in", res = 500)
  par(mfrow = c(3, 3), mar = c(0, 0, 0, 0))
  while (temp_n < 150) {
    
    ## Make grid based on latlon coords
    temp_grid <- sp::makegrid(x = grid_pts_aea, cellsize = temp_res)
    temp_grid <- sp::SpatialPoints(coords = coordinates(temp_grid),
                                   proj4string = aea_crs)
    temp_grid <- crop(x = temp_grid, y = chukchi_mask)
    temp_n <- length(temp_grid)
    
    # plot(chukchi_mask); points(temp_grid, pch = 16)
    
    grid_idx <- apply(X = gDistance(spgeom1 = grid_pts_aea,
                                    spgeom2 = temp_grid,
                                    byid = TRUE),
                      MARGIN = 1,
                      FUN = which.min)
    
    plot(chukchi_mask, asp = 1)
    box()
    legend("bottomright", legend = paste0("n = ", temp_n), bty = "n", cex = 2)
    points(temp_grid,  pch = 16, cex = 1)
    
    sys_settings <- rbind(sys_settings,
                          data.frame(n = temp_n,
                                     res = temp_res))
    
    temp_res <- temp_res - 5000
    sys_idx <- c(sys_idx, list(grid_idx))
  }
  
  # dev.off()
}

##################################
## Choose a gear to simulate surveys
##################################
igear = c("otter", "beam")[1]
true_index <- get(paste0("true_index_", igear))
true_index <- true_index[nrow(true_index), ]

##################################
## Synthesize simulated densities across species for a given igear.
## For now, only simulate the most recent year of data for a given igear.
##################################
n_spp <-  get(paste0("n_spp_", igear))
spp_list <-  get(paste0("spp_list_", igear))
iyear <- c("beam" = "2019", "otter" = "2012")[igear]
ms_dens <- array(dim = c(n_spp, n_cells, 1000), 
                 dimnames = list(spp_list, NULL))

for (ispp in 1:n_spp) { ## Loop over species -- start
  load(paste0("results/chukchi_", igear,
              "/vast_fits/", spp_list[ispp], "/simulated_densities/sim_data_",
              iyear, "_simtype1.RData"))
  
  ms_dens[ispp, , ] <- get(paste0("sim_data_", iyear, "_simtype1"))
} ## Loop over species -- end

##################################
## Simulate Simple Random Designs at varying sampling efforts
##################################
target_n <- seq(from = 50, by = 10, to = 170)
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

##################################
## Simulate Systematic Random Designs at varying sampling efforts
## Variance is calculated using a local calculation
##################################
index_sys <- cv_sys <- rb_index_sys <-
  array(dim = c(nrow(sys_settings), n_spp, n_iters),
        dimnames = list(paste0("sys_n = ", sys_settings$n), spp_list, NULL))

set.seed(234)
for (isample in 1:nrow(sys_settings)) { ## Loop over survey effort -- start
  
  ## Calculate the nearest neighbor stations for each station
  distances <- as.matrix(dist(grid_pts_aea[sys_idx[[isample]], ]@coords))
  min_distances <- apply(X = distances, 
                         MARGIN = 1, 
                         FUN = function(x) {
                           y = round(x[x > 0])
                           return(names(which( y == min(y) )))
                         }) 
  
  for (iter in 1:n_iters) { ## Loop over iterations -- start
    
    neighbor_var <- matrix(nrow = length(min_distances), ncol = n_spp,
                           dimnames = list(NULL, spp_list))
    
    for (icell in 1:length(min_distances) ){  ## Loop over stations -- start
      
      ## Pull the station indices of the station and its neighbors
      idx <- as.integer(c(icell, min_distances[[icell]]))
      
      ## Pull the location of the idx on the grid from sys_idx and get the 
      ## simulated densities
      neighbors <- ms_dens[, sys_idx[[isample]][idx], iter]
      
      ## Calculate the neighborhood variance of the station
      neighbor_var[icell, ] <-  apply(X = neighbors, 
                                      MARGIN = 1, 
                                      FUN = function(x) mean((x - mean(x))^2))
    } ## Loop over stations -- end
    
    ## Calculate survey statistics
    sys_sd_LO5 <- sqrt(colSums(neighbor_var) / length(min_distances)^2)
    sys_mean <- rowMeans(ms_dens[, sys_idx[[isample]], iter])
    
    ## Record survey statistics
    index_sys[isample, , iter] <- sys_mean * total_area
    cv_sys[isample, , iter] <- sys_sd_LO5 / sys_mean
    rb_index_sys[isample, , iter] <- 
      100 * (index_sys[isample, , iter] - true_index) / true_index
    
  }  ## Loop over iterations -- end
  print(paste0("Finished with SYS n = ", sys_settings$n[isample]))
} ## loop over survey effort -- end

##################################
## Multi Species STRS Optimization
##################################
index_ms_strs <- cv_ms_strs <- rb_index_ms_strs <-
  array(dim = c(length(target_n), n_spp, n_iters),
        dimnames = list(paste0("ms_strs_n = ", target_n), spp_list, NULL))
n_strata <- 5

load(paste0("results/chukchi_", igear,
            "/survey_opt/Str_", n_strata, "/result_list.RData"))
load(paste0("results/chukchi_", igear,
            "/survey_opt/Str_", n_strata, "/allocations.RData"))

set.seed(231)
for (itarget in target_n){ ## Loop over sampling effort -- start
  
  ## STRS settings
  solution <- result_list$sol_by_cell
  allocation <- as.integer(subset(ms_sample_allocations, 
                                  n == itarget)[, paste0("Str_ ", 1:n_strata)])
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
} ## Loop over sampling effort -- end

##################################################
####   Calculate Performance Metrics
##################################################  
for (idesign in c("srs", "sys", "ms_strs")) { ## Loop over design type -- start
  
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

##################################################
####   Save
##################################################  
save(file = paste0("results/chukchi_", igear, "/survey_sim_results.RData"),
     list = c(as.vector(sapply(X = c("true_cv_", "rrmse_cv_", 
                                     "index_", "cv_", 
                                     "rb_index_", "rb_cv_"), 
                               FUN = function(x) 
                                 paste0(x, c("srs", "sys", "ms_strs")))),
              "sys_settings", "sys_idx", "target_n"))
