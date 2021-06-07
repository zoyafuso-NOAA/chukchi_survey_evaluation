###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

##################################################
####   Set up directory
##################################################
wd <- "C:/Users/zack.oyafuso/Desktop/Arctic/"

##################################################
####   Load packages
##################################################
library(raster)
library(sp)
library(rgdal)

##################################################
####   Import Data
##################################################
chukchi_grid <- read.csv(paste0(wd, 
                                "spatial_data/BS_Chukchi_extrapolation_grids/",
                                "ChukchiThorsonGrid.csv"))

##################################################
####  Import land data
##################################################
AK <- rgdal::readOGR(paste0(wd, "spatial_data/land_shapefiles/AKland.shp"))
# AK <- sp::spTransform(AK, CRSobj = crs("+proj=longlat") )

##################################################
####   Covariate names from J. Marsh
##################################################
covariate_names <- unique(tools::file_path_sans_ext(
  dir(paste0(wd, "covariate_data/"))))[-1]

##################################################
####   Import covariates and turn into a raster brick
##################################################
covariate_brick <- brick()

for (icovar in 1:length(covariate_names)) {
  covariate_brick = brick(
    x = list(covariate_brick,
             raster(paste0(wd, "covariate_data/", 
                           covariate_names[icovar], ".grd")) ))
}
names(covariate_brick) <- covariate_names
plot(covariate_brick)

##################################################
####   Import bathymetry covarites from J. Purtle
##################################################
covariate_names_bathy <- c("depth" = "arctic1km_new",
                           "bpi_1km" = "bpi65_1km",
                           "slope" = "slope_1km")

# covariate_bathy_brick <- brick()

for (icovar in 1:length(covariate_names_bathy)) {
  covariate_brick = raster::brick(
    x = list(covariate_brick,
             raster::raster(paste0(wd, 
                                   "covariate_data/Arctic Bathymetry 1 km_NEW/", 
                                   covariate_names_bathy[icovar], "/")) 
    ))
}

idx <- (1 + length(covariate_names)):((1+length(covariate_names)) + 
                                        length(covariate_names_bathy) - 1)

names(covariate_brick)[idx] <- names(covariate_names_bathy)

par(mfrow = c(4, 5), mar = c(1, 1, 1, 4))
for(icovar in c(covariate_names, names(covariate_names_bathy)) ) {
  plot(covariate_brick[[icovar]], axes = F)
  mtext(text = icovar, side = 3, line = -1.5)
}


##################################################
####   COnvert chukchi grid cells to a spatial dataframe with aea projection
##################################################
chukchi_pts_aea <- sp::spTransform(
  x = SpatialPointsDataFrame(coords = chukchi_grid[, c("Lon", "Lat")],
                             proj4string = CRS("+proj=longlat"),
                             data = chukchi_grid), 
  CRSobj = covariate_brick@crs)

##################################################
####   Extract covariate data from raster brick to each grid cell
##################################################
chukchi_pts_aea@data[, c("E_km", "N_km")] <- chukchi_pts_aea@coords
chukchi_pts_aea@data[, c(covariate_names, names(covariate_names_bathy))] <- 
  raster::extract(x = covariate_brick, 
                  y = chukchi_pts_aea)

summary(chukchi_pts_aea@data)

##################################################
####   For the missing non-overlapping grid cells, choose the closest 
####   raster cell with covariate data
##################################################

missing_idx_matrix <- matrix(data = FALSE,
                             nrow = nrow(chukchi_pts_aea),
                             length(names(covariate_brick)),
                             dimnames = list(NULL, names(covariate_brick)))

for (icovar in names(covariate_brick)) {
  r <- covariate_brick[[icovar]] #example raster set
  missing_idx <- is.na(chukchi_pts_aea@data[, icovar])
  sampled = apply(X = chukchi_pts_aea@data[missing_idx, c("E_km", "N_km")] , 
                  MARGIN = 1, 
                  FUN = function(xy) 
                    which.min(replace(distanceFromPoints(r, xy), 
                                      is.na(r), 
                                      NA)) )
  
  if (length(sampled) > 0) {
    chukchi_pts_aea@data[missing_idx, icovar] <- r[sampled]
    missing_idx_matrix[missing_idx, icovar] <- TRUE
  }
  
  print(paste("Done with", icovar))
}

summary(chukchi_pts_aea)

##################################################
####   Plot chukchi covariates with missing data
##################################################

par(mfrow = c(3, 7), mar = c(1, 1, 1.5, 5))
for(icovar in c(covariate_names, names(covariate_names_bathy))) {
  chukchi_ras <- raster::raster(x = chukchi_pts_aea, 
                                resolution = 4000)
  chukchi_ras <- raster::rasterize(x = chukchi_pts_aea, 
                                   y = chukchi_ras, 
                                   field = icovar)
  plot(chukchi_ras, main = icovar, axes = F)
  points(chukchi_pts_aea@data[missing_idx_matrix[, icovar], c("E_km", "N_km")],
         pch = 16, cex = 0.1)
  # plot(AK, add = T, col = "tan", border = FALSE)
}

round(cor(chukchi_pts_aea@data[, names(covariate_brick)], 
          use = "complete.obs"), 
      digits = 2)

##################################################
####   Save
##################################################
save(list = c("chukchi_pts_aea", "covariate_brick"),
     file = paste0(wd, "covariate_data/covariate_data.RData"))

