###############################################################################
## Project:      Chukchi STRS Design Solutions
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Show the 100 station solution for the Chukchi otter trawl
##                   STRS design optimization with 3, 4, and 5 strata
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)

##################################################
####   Constants and result objects
##################################################
load("data/survey_opt_data/optimization_data.RData")
akland <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")

strata <- 3:5
solution <- matrix(nrow = n_cells, ncol = 3, 
                   dimnames = list(NULL, paste0("Str_", strata)))
allocations <- list()

##################################################
####   Collect solution for each level of strata
##################################################
for (istrata in strata) {
  result_dir <- paste0("results/otter_trawl/survey_opt - ", 
                       istrata, " STRATA/MS/")
  nruns <- length(dir(result_dir))
  load(paste0(result_dir, "Run_", nruns, "/result_list.RData"))
  
  solution[, paste0("Str_", istrata)] <- result_list$sol_by_cell
  allocations[[paste0("Str_", istrata)]] <- result_list$sum_stats$Allocation
}

##################################################
####   Plot each solution
##################################################
{
  png(filename = "presentations/results_12_07_2021/strata_sols.png",
     width = 8, height = 4, units = "in", res = 500)
  par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))
  for (istrata in strata) {
    chukchi_shp <- sp::SpatialPointsDataFrame(coords = grid_pts@coords, 
                                              proj4string = aea_crs, 
                                              data = data.frame(solution))
    
    chukchi_ras <- raster::raster(chukchi_shp, res = 3700)
    chukchi_ras <- raster::rasterize(x = chukchi_shp, 
                                     y = chukchi_ras, 
                                     field = paste0("Str_", istrata))
    
    sample_vec <- c()
    for(ilev in 1:istrata) {
      sample_vec <- 
        c(sample_vec,
          sample(x = which(solution[, paste0("Str_", istrata)] == ilev),
                 size = allocations[[paste0("Str_", istrata)]][ilev]) )
    }
    
    samp_pts <- sp::SpatialPoints( coords = grid_pts@coords[sample_vec, ], 
                                   proj4string = aea_crs)
    
    image(chukchi_ras, 
          col = rev(RColorBrewer::brewer.pal(name = "Paired", n = istrata)),
          axes = F, ann = F, asp = 1)
    plot(akland, add = TRUE, col = "tan", border = FALSE)
    points(samp_pts,
           pch = 16, cex = 0.5)
    legend("topright", legend = paste(istrata, "Strata"),
           bty = "n", cex = 1.5, text.font = 2)
    box()
  }
  dev.off()
}
