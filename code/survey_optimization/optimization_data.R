###############################################################################
## Project:      Chukchi optimization 
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Constants used in survey optimization
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(sp)
library(rgdal)
library(rgeos)
library(raster)

##################################################
#### Import Chukchi Grid and AK land   
##################################################
chukchi_grid <- read.csv(paste0("data/spatial_data/",
                                "BS_Chukchi_extrapolation_grids/",
                                "ChukchiThorsonGrid.csv"))
AKland <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")

aea_crs <- raster::crs(AKland)
latlon_crs <- sp::CRS("+proj=longlat +datum=WGS84")
n_cells <- nrow(chukchi_grid)

##################################################
#### Calculate minimum distance to land for each grid point
##################################################
grid_pts <- sp::SpatialPoints(coords = chukchi_grid[, c("Lon", "Lat")], 
                              proj4string = latlon_crs)
grid_pts <- sp::spTransform(x = grid_pts, CRSobj = aea_crs)

dist_to_shore <- apply(X = rgeos::gDistance(spgeom1 = grid_pts, 
                                            spgeom2 = AKland, 
                                            byid = TRUE), 
                       MARGIN = 2, 
                       FUN = min)

xrange <- bbox(grid_pts)[1, ]
yrange <- bbox(grid_pts)[2, ]

##################################################
#### Species included
##################################################
spp_list <- c("Alaska plaice", "Arctic cod", "Bering flounder", "saffron cod",
              "snow crab", "yellowfin sole")
n_spp <- length(spp_list)
names(spp_list) <- paste0("Y", 1:n_spp)

##################################################
#### Species optimization data input
#### X1 and X2 are stratum variables latitude (X1) and distance from shore (X2)
#### Y1, Y2, ..., Yn_spp are density predictions from VAST and 
#### Y1_SUM_SQ, ..., is the square of the density prediction
##################################################
frame_df <- data.frame(domainvalue = 1,
                       id = 1:nrow(chukchi_grid),
                       WEIGHT = 1,
                       X1 = chukchi_grid$Lat,
                       X2 = dist_to_shore)

true_index <- c()

for (ispp in 1:n_spp) {
  load(paste0("results/otter_trawl/", spp_list[ispp], "/fit.RData"))
  frame_df[, paste0("Y", ispp)] <- fit$Report$D_gct[, 1, 2]
  frame_df[, paste0("Y", ispp, "_SQ_SUM")] <- fit$Report$D_gct[, 1, 2]^2
  
  true_index[ispp] <- fit$Report$Index_ctl[1, 2, 1]
}

names(true_index) <- paste0("Y", 1:n_spp)

##################################################
#### Save
##################################################
save(list = c("aea_crs", "latlon_crs", "xrange", "yrange", "grid_pts", 
              "n_cells", "spp_list", "n_spp", "frame_df", "true_index"),
     file = "data/survey_opt_data/optimization_data_otter.RData")

