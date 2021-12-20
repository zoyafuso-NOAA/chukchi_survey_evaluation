###############################################################################
## Project:      Determine 2012 BTS Chukchi Stations within 3 nautical miles
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Calculate closest distance of each station to the AK polygon
##                        Query which stations are within 3 nm (5556 meters)
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(rgdal)
library(sp)
library(raster)

boundary <- 3 * 1852 ## Conversion from 3 nm to km

##################################################
####   Import Datsets, Crop AKland to just those in the Chukchi
##################################################
AKland <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
aea_crs <- raster::crs(AKland)
latlon_crs <- sp::CRS("+proj=longlat +datum=WGS84")

AK_cropped <- raster::crop(x = AKland, 
                           y = extent(-684558.5, -664.9984, 
                                      1779066,  2770210
                                      ))

AK_buffer <- raster::buffer(x = AK_cropped, width = boundary)

##################################################
####   Subset 2012 Chukchi otter trawldata
##################################################
bts <- subset(x = read.csv(paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                  "AK_BTS_Arctic_processed_wide.csv")),
              subset = SURVEY_NAME == "2012 Chukchi Sea Trawl Survey" &
                GEAR_CAT == "otter")

##################################################
####   Calculate closest distance from each station location to land
##################################################
grid_pts <- sp::SpatialPoints(coords = bts[, c("MEAN_LONGITUDE", 
                                               "MEAN_LATITUDE")], 
                              proj4string = latlon_crs)
grid_pts <- sp::spTransform(x = grid_pts, CRSobj = aea_crs)

dist_to_shore <- apply(X = rgeos::gDistance(spgeom1 = grid_pts, 
                                            spgeom2 = AKland, 
                                            byid = TRUE), 
                       MARGIN = 2, 
                       FUN = min)

##################################################
####   Calculate closest distance from each station location to land
##################################################
par(mar = c(0,0,0,0))
plot(AK_buffer, col = "red", border = F)
plot(AK_cropped, col = "tan", border = F, add = TRUE); 
points(grid_pts, pch = 16, cex = 0.5)
