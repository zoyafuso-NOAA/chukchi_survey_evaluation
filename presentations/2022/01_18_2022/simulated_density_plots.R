###############################################################################
## Project:       Simulated Density Plots
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Sim type 1 and 3
###############################################################################
rm(list = ls())

##################################################
#### Import relevant libraries   
##################################################
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)
library(gap)
library(plotrix)

##################################################
#### Import Data and Constants
##################################################
load("data/survey_opt_data/optimization_data.RData")
AK_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
cropped_extent <- extent(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- raster::crop(x = AK_land, y = cropped_extent)

##################################################
#### Synthesize density predictions for each gear and species
##################################################
for (igear in c("otter", "beam")) {
  
  ## Years that gear was used
  years <- c("otter" = 2012, "beam" = 2019)[igear]
  
  ## Species that gear sampled
  species <- list.dirs(paste0("results/chukchi_", igear, "/vast_fits/"), 
                       full.names = FALSE, 
                       recursive = FALSE)
  
  ## Temporary result object
  ms_dens <- array(data = NA,
                   dim = c(n_cells, length(species), 10),
                   dimnames = list(NULL, species, NULL))
  
  for (ispp in species) {
    ## Load VAST results 
    load(paste0("results/chukchi_", igear, 
                "/vast_fits/", ispp, "/sim_data_", years, 
                "_simtype3_iter1to500.RData"))
    
    ## Append to temporary result object
    ms_dens[, ispp, ] <- get(paste0("sim_data_", years, "_simtype3"))[, 1:10]
  }
  
  ## attach gear type with ms_dens variable name
  assign(x = paste0("ms_dens_", igear), value = ms_dens)
}

##################################################
#### Plot
##################################################
col_otter <- brewer.pal(n = 9, name = "Blues")
col_beam <- brewer.pal(n = 9, name = "Greens")

for (ispp in c("Alaska plaice", "Arctic cod", "Bering flounder", 
               "saffron cod", "snow crab", "yellowfin sole")){
  
  ## Setup layout of 
  par(mar = c(0, 0, 0, 0))
  
  available_gears <- switch(ispp,
                            "Alaska plaice" =  "otter", 
                            c("otter", "beam"))
  
  for (igear in available_gears) {
    
    png(filename = paste0("presentations/2022/01_18_2022/sim3_density_", 
                          igear, "_", ispp, ".png"),
        width = 100, units = "mm", res = 500, family = "serif",
        height = 110)
    par(mar = c(0, 0, 0, 0), mfrow = c(3, 3))
    ## Years that gear was used
    years <- c("otter" = 2012, "beam" = 2019)[igear]
    
    for (iter in 1:9) {
      temp_density <- get(paste0("ms_dens_", igear))[, ispp, iter]
      density_cuts <- round(c(0, 1, 
                              quantile(x = temp_density[temp_density > 1],
                                       probs = c(0.25, 0.5, 0.75)),
                              max(temp_density) + 1), 1)
      assign(x = paste0("density_cuts_", igear), value = density_cuts)
      
      plot_this <- data.frame(X1 = temp_density)
      
      plot_this$X1 <- as.numeric(base::cut(x = plot_this$X1, 
                                           breaks = density_cuts, 
                                           labels = 1:(length(density_cuts) - 1),
                                           include.lowest = TRUE,
                                           right = FALSE))
      chukchi_shp <- sp::SpatialPointsDataFrame(coords = grid_pts@coords, 
                                                proj4string = aea_crs,
                                                data = plot_this)
      chukchi_ras <- raster::raster(x = chukchi_shp, resolution = 3700)
      chukchi_ras <- raster::rasterize(x = chukchi_shp,
                                       y = chukchi_ras, 
                                       field = "X1")
      
      image(chukchi_ras, asp = 1, col = get(paste0("col_", igear)), 
            axes = F, ann = F)
      box()
      
      ## Plot land
      plot(ak_land_cropped, add = T, border = F, col = "tan")
    
      ## Legend
      plotrix::color.legend(
        xl = extent(grid_pts)[1] + diff(extent(grid_pts)[1:2]) * 0.8, 
        xr = extent(grid_pts)[1] + diff(extent(grid_pts)[1:2]) * 0.9,
        yb = extent(grid_pts)[3] + diff(extent(grid_pts)[3:4]) * 0.025, 
        yt = extent(grid_pts)[3] + diff(extent(grid_pts)[3:4]) * 0.60, 
        legend = c("<1", round(get(paste0("density_cuts_", igear)))[-(1:2)]),
        rect.col = get(paste0("col_", igear)),
        gradient = "y", 
        align = "rb",
        cex = 0.4)
    }
    
    dev.off()
  }
  
}
