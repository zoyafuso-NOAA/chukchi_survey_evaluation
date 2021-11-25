###############################################################################
## Project:       For 23 Nov 2021 Updates
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Otter Trawl STRS optimization preliminary results
###############################################################################
rm(list = ls())

##################################
##  Import Library
##################################
library(raster)
library(sp)

##################################
## Load Data 
##################################
load("results/otter_trawl/survey_opt/survey_res.RData")
load("results/otter_trawl/survey_opt/SS/ss_knit_res.RData")
load("data/survey_opt_data/optimization_data_otter.RData")
load("results/otter_trawl/survey_opt/MS/Run_14/result_list.RData")

##################################
## Performance Metrics
##################################
png(filename = "presentations/results_11_23_2021/rb.png",
    width = 6, height = 6, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(3, 3, 1.5, 1), oma = c(2, 2, 0, 0))
for (ispp in 1:n_spp) {
  boxplot(100 * rb[ispp, , ], ylim = c(-100, 200), las = 1, main = spp_list[ispp])
  abline(h = 0)
}
mtext(side = 1, text = "Survey Design (n = 71)", outer = T, font = 2)
mtext(side = 2, text = "Percent Bias", outer = T, font = 2)
dev.off()

png(filename = "presentations/results_11_23_2021/cvs.png",
    width = 6, height = 6, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(3, 3, 1.5, 1), oma = c(2, 2, 0, 0))
for (ispp in 1:n_spp) {
  boxplot(cv[ispp, , ], 
          ylim = c(0, 
                   min(1, 
                       max(c(max(true_cv[ispp, ]), 
                             max(cv[ispp, , ]))) )), 
          las = 1, main = spp_list[ispp])
  points(x = 1:4, 
         y = ifelse(test = true_cv[ispp, ] > 1, 
                    yes = 0.9, 
                    no = true_cv[ispp, ]), 
         pch = "*", cex = 3, 
         col = ifelse(test = true_cv[ispp, ] > 1, 
                      yes = "blue", 
                      no = "red") )
}

legend("topright", legend = "True CV", col = "red", pch = "*", pt.cex = 3)
mtext(side = 1, text = "Survey Design (n = 71)", outer = T, font = 2)
mtext(side = 2, text = "Survey CVs", outer = T, font = 2)
dev.off()

##################################
## Simulated Densities
##################################
ms_dens <- array(dim = c(n_spp, n_cells, 500), 
                 dimnames = list(spp_list, NULL, NULL))
for (ispp in 1:n_spp) {
  ## Load simulated densities
  load(paste0("results/otter_trawl/", spp_list[ispp], "/simulated_data.RData"))
  ms_dens[ispp, , ] <- sim_data
}

png(filename = "presentations/results_11_23_2021/sim_data_vis.png",
    width = 4, height = 6, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
for (ispp in 1:n_spp) {
  for (iter in 100){
    goa <- sp::SpatialPointsDataFrame(
      coords = grid_pts@coords,
      data = data.frame(Str_no = ms_dens[ispp, , iter]), 
      proj4string = crs(grid_pts) )
    goa_ras <- raster::raster(x = goa,
                              resolution = 5000)
    goa_ras <- raster::rasterize(x = goa,
                                 y = goa_ras,
                                 field = "Str_no")
    plot(goa_ras, axes = F)
    legend("bottomright", legend = spp_list[ispp], bty = "n")
  }
}
dev.off()

##################################
## Single Species Optimizations
##################################
png(filename = "presentations/results_11_23_2021/ss_solns.png",
    width = 4, height = 6, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
for (ispp in 1:n_spp) {
  goa <- sp::SpatialPointsDataFrame(
    coords = grid_pts@coords,
    data = data.frame(Str_no = ss_solutions[, ispp ]), 
    proj4string = crs(grid_pts) )
  goa_ras <- raster::raster(x = goa,
                            resolution = 5000)
  goa_ras <- raster::rasterize(x = goa,
                               y = goa_ras,
                               field = "Str_no")
  plot(goa_ras, axes = F)
  
  sample_vec <- c()
  for(istrata in 1:nrow(ss_allocations)) {
    sample_vec <- c(sample_vec,
                    sample(x = which(ss_solutions[, ispp] == istrata),
                           size = ss_allocations[istrata, ispp]) )
  }
  
  samp_pts <- sp::SpatialPoints(
    coords = grid_pts@coords[sample_vec, ], 
    proj4string = crs(grid_pts))
  
  points(samp_pts,
         pch = 16, cex = 0.5)
  
  legend("bottomright", legend = spp_list[ispp], bty = "n")
}
dev.off()

##################################
## Multispecies Solution
##################################
png(filename = "presentations/results_11_23_2021/ms_solns.png",
    width = 6, height = 6, units = "in", res = 500)
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))

  goa <- sp::SpatialPointsDataFrame(
    coords = grid_pts@coords,
    data = data.frame(Str_no = result_list$sol_by_cell), 
    proj4string = crs(grid_pts) )
  goa_ras <- raster::raster(x = goa,
                            resolution = 5000)
  goa_ras <- raster::rasterize(x = goa,
                               y = goa_ras,
                               field = "Str_no")
  plot(goa_ras, axes = F)
  
  sample_vec <- c()
  for(istrata in 1:nrow(ss_allocations)) {
    sample_vec <- c(sample_vec,
                    sample(x = which(result_list$sol_by_cell == istrata),
                           size = result_list$sum_stats$Allocation) )
  }
  
  samp_pts <- sp::SpatialPoints(
    coords = grid_pts@coords[sample_vec, ], 
    proj4string = crs(grid_pts))
  
  points(samp_pts, pch = 16, cex = 0.5)

dev.off()
