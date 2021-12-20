###############################################################################
## Project:       Extract Static and Dynamic Covariates
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Extract covariate data to the trawl locations 
##                        and interpolation grid
##
##                Static Variables: Depth, SLope and BPI at 1 km scale
##                Dynamic Variables: bottom temperature, salinity, 
##                        current direction (horizontal and vertical direction)
###############################################################################
rm(list = ls())

##################################################
####   Load packages
##################################################
library(raster)
# library(raster)
# library(rgdal)
# library(ncdf4)

##################################################
####   Import Interpoloation Grid
##################################################
chukchi_grid <- 
  read.csv(paste0("data/spatial_data/BS_Chukchi_extrapolation_grids/",
                  "ChukchiThorsonGrid.csv"))

##################################################
####   Import Static and Dynamic Covariates
##################################################
for (ifile in  c("static_covariate_brick", 
                 "temp_brick", "salt_brick", "u_brick", "v_brick")) {
  assign(x = ifile, 
         value = brick(x = paste0("data/covariate_data/", ifile, ".grd")))
}


##################################################
####   Import Otter Trawl Data
##################################################
otter <- read.csv(paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                         "AK_BTS_Arctic_processed_wide.csv"))

##################################################
####   Convert chukchi grid cells (chukchi_pts_aea) and the station locations
####   of the otter trawl (otter_pts_aea) to a spatial dataframe with aea 
####   projection
##################################################
chukchi_pts_aea <- sp::spTransform(
  x = sp::SpatialPointsDataFrame(coords = chukchi_grid[, c("Lon", "Lat")],
                                 proj4string = CRS("+proj=longlat"),
                                 data = chukchi_grid),
  CRSobj = crs(temp_brick))

otter_pts_aea <- sp::spTransform(
  x = sp::SpatialPointsDataFrame(coords = otter[, c("MEAN_LONGITUDE", 
                                                    "MEAN_LATITUDE")],
                                 proj4string = CRS("+proj=longlat"),
                                 data = otter),
  CRSobj = crs(temp_brick))

##################################################
####   Extract static covariate data from raster brick to each grid cell
####   in chukchi_pts_aea dn otter_pts_aea
##################################################
chukchi_pts_aea@data[, c("E_km", "N_km")] <- chukchi_pts_aea@coords
chukchi_pts_aea@data[, names(static_covariate_brick)] <-
  raster::extract(x = static_covariate_brick,
                  y = chukchi_pts_aea)

otter_pts_aea@data[, c("E_km", "N_km")] <- otter_pts_aea@coords
otter_pts_aea@data[, names(static_covariate_brick)] <-
  raster::extract(x = static_covariate_brick,
                  y = otter_pts_aea)

##################################################
####   Extract dynamic covariate data from raster brick to each grid cell
####   in chukchi_pts_aea dn otter_pts_aea
##################################################
for (ivar in c("temp", "u", "v", "salt")) {
  for (iyear in c(1990, 2012)) {
    temp_ras <- get(paste0(ivar, "_brick"))[[paste0("X", iyear, "_Aug_Avg")]]
    chukchi_pts_aea@data[, paste0(ivar, "_", iyear)] <-
      raster::extract(x = temp_ras,
                      y = chukchi_pts_aea)
  }
  
  for (idate in unique(otter_pts_aea@data$DATE)) {
    temp_date <- gsub(x = paste0("X", idate), pattern = "-", replacement = ".")
    temp_ras <- get(paste0(ivar, "_brick"))[[temp_date]]
    idx <- otter_pts_aea@data$DATE == idate
    
    otter_pts_aea@data[idx, paste0(ivar)] <-
      raster::extract(x = temp_ras,
                      y = otter_pts_aea[idx, ]) 
  }
}

##################################################
####   For the missing non-overlapping grid cells, choose the closest
####   raster cell with static covariate data
##################################################
missing_idx_matrix <- matrix(data = FALSE,
                             nrow = nrow(chukchi_pts_aea),
                             length(names(static_covariate_brick)),
                             dimnames = list(NULL, names(static_covariate_brick)))

for (icovar in names(static_covariate_brick)[1]) {
  r <- static_covariate_brick[[icovar]] #example raster set
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

##################################################
####   For the missing non-overlapping grid cells, choose the closest
####   raster cell with static covariate data
##################################################
missing_idx_matrix <- 
  matrix(data = FALSE,
         nrow = nrow(chukchi_pts_aea),
         ncol = length(paste0(rep(c("temp", "u", "v", "salt"), 
                                  each = 2), 
                              "_", 
                              rep(c(1990, 2012), 
                                  times = 4))),
         dimnames = list(NULL, paste0(rep(c("temp", "u", "v", "salt"), 
                                          each = 2), 
                                      "_", 
                                      rep(c(1990, 2012), 
                                          times = 4))))

for (ivar in c("temp", "u", "v", "salt")) {
  for (iyear in c(1990, 2012)) {
    r <- get(paste0(ivar, "_brick"))[[paste0("X", iyear, "_Aug_Avg")]] #temp raster
    missing_idx <- is.na(chukchi_pts_aea@data[, paste0(ivar, "_", iyear)])
    sampled = apply(X = chukchi_pts_aea@data[missing_idx, c("E_km", "N_km")] ,
                    MARGIN = 1,
                    FUN = function(xy)
                      which.min(replace(distanceFromPoints(r, xy),
                                        is.na(r),
                                        NA)) )
    
    if (length(sampled) > 0) {
      chukchi_pts_aea@data[missing_idx, paste0(ivar, "_", iyear)] <- r[sampled]
      missing_idx_matrix[missing_idx, paste0(ivar, "_", iyear)] <- TRUE
    }
    
    print(paste0("X", iyear, "_Aug_Avg"))
  }
}

##################################################
####   Combine the station and grid data
##################################################
wide_covar_df <- cbind(type = "grid", 
                       chukchi_pts_aea@data[, c("Lon", "Lat", "depth", 
                                                paste0(rep(c("temp", "salt", 
                                                             "u", "v"), 
                                                           each = 2), 
                                                       "_", 
                                                       rep(c(1990, 2012), 
                                                           times = 4)))])

temp_otter <- cbind(type = "station",
                    otter_pts_aea@data[, c("MEAN_LONGITUDE", "MEAN_LATITUDE",
                                           "YEAR",
                                           "depth", "temp", "salt", "u", "v")])
names(temp_otter)[2:4] <- c("Lon", "Lat", "year")

covar_df <- data.frame()
for (iyear in c(1990, 2012)) {
  temp_df <- cbind(year = iyear,
                   wide_covar_df[, c("type", "Lon", "Lat", "depth",
                                     paste0(c("temp", "salt", "u", "v"), 
                                            "_", iyear))])
  names(temp_df)[6:9] <- c("temp", "salt", "u", "v")
  covar_df <- rbind(covar_df, temp_df)
}
covar_df <- rbind(covar_df, temp_otter[, names(covar_df)])

##################################################
####   Save
##################################################
save(file = "data/covariate_data/covariate_df.RData",
     list = "covar_df")
