###############################################################################
## Project:       Wrangle Static and Dynamic Chukchi Sea Covariates
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Static: Depth, SLope and BPI at 1 km scale
##                Dynamic: bottom temperature, salinity, current direction 
##                         (horizontal and vertical direction )
##                Done in R Version 4.0.3
###############################################################################
rm(list = ls())

##################################################
####   Load packages
##################################################
library(raster)
library(sp)
library(rgdal)
library(ncdf4)

##################################################
####   Location of dynamic variables, Lewis Barnett's G drive
##################################################
dyn_covar_dir <- "G:/Barnett/PAROMS/Arctic4/averages/"

##################################################
####   Import grid of the dynamic covariates 
####   (Note: different from the interpolation grid)
##################################################
arctic_grid <- ncdf4::nc_open(paste0(dirname(dyn_covar_dir), 
                                     "/grid_Arctic_4.nc"))

arctic_lon <- ncvar_get(arctic_grid, "lon_rho")
arctic_lon <- ifelse(arctic_lon > 180, -360 + arctic_lon, arctic_lon)
arctic_lat <- ncvar_get(arctic_grid, "lat_rho")

##################################################
####   Import locations of otter trawls
##################################################
otter_df <- read.csv(file = paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                   "AK_BTS_Arctic_processed_wide.csv"))

##################################################
####  Import land data
##################################################
AK_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")

##################################################
####   Import bathymetry covariates from J. Purtle and put into a RasterBrick
##################################################
covariate_names_bathy <- c("depth" = "arctic1km_new",
                           "bpi_1km" = "bpi65_1km",
                           "slope" = "slope_1km")

static_covariate_brick <- raster::brick()

for (icovar in 1:length(covariate_names_bathy)) {
  static_covariate_brick = raster::brick(
    x = list(static_covariate_brick,
             raster::raster(paste0("data/covariate_data/",
                                   "Arctic Bathymetry 1 km_NEW/", 
                                   covariate_names_bathy[icovar], "/")) 
    ))
}

names(static_covariate_brick) <- names(covariate_names_bathy)

##################################################
####   Import dynamic covariates
##################################################
temp_brick <- salt_brick <- u_brick <- v_brick <- NULL

unique_dates <- unique(otter_df$DATE)
lon_idx <- 310:600; lat_idx <- 51:550; depth_layer <- 1

for (idate in unique_dates){
  iyear <- substr(x = idate, start = 1, stop = 4)
  imonth <- substr(x = idate, start = 6, stop = 7)
  iday <- substr(x = idate, start = 9, stop = 10)
  
  ifile <- paste0(dyn_covar_dir, 
                  "arctic4_avg_", iyear, "-", imonth, "-", iday, "T.nc")
  
  temp_nc <- ncdf4::nc_open(filename = ifile)
  
  for (ivar in c("temp", "salt", "u", "v")) {
    temp_var <- ncvar_get(temp_nc, ivar)
    
    goa <- sp::SpatialPointsDataFrame(
      coords = cbind(lon = as.numeric(arctic_lon[lon_idx, lat_idx]),
                     lat = as.numeric(arctic_lat[lon_idx, lat_idx])),
      data = data.frame(var = as.numeric(temp_var[lon_idx, 
                                                  lat_idx, 
                                                  depth_layer])))
    
    goa_ras <- raster::raster(x = goa,
                              nrows = diff(range(lon_idx)),
                              ncols = diff(range(lat_idx)))
    goa_ras <- raster::rasterize(x = goa,
                                 y = goa_ras,
                                 field = "var",
                                 na.rm = FALSE)
    
    projection(goa_ras) <- CRS("+proj=longlat +datum=WGS84")
    goa_ras <- projectRaster(goa_ras, crs = crs(AK_land))
    
    if (idate == unique_dates[1]) {
      assign(value = raster::brick(goa_ras), 
             x = paste0(ivar, "_brick"))
    }
    if (idate != unique_dates[1]) {
      assign(value = raster::brick(list(get(paste0(ivar, "_brick")),
                                        goa_ras)), 
             x = paste0(ivar, "_brick"))
    }
  }
  
  print(paste0("Finished with ", idate))
}

names(temp_brick) <- names(v_brick) <- names(u_brick) <- names(salt_brick) <- 
  unique_dates

##################################################
####   Attach Average August values for use in the extrapolation grids
##################################################
years_to_extract <- c(1990, 2012)
month_to_extract <- 8 #August

for (iyear in c(1990, 2012)) {
  ifile <- paste0(dirname(dyn_covar_dir), "/months/",
                  "arctic4_", iyear, "_0", month_to_extract, ".nc")
  
  temp_nc <- ncdf4::nc_open(filename = ifile)
  
  for (ivar in c("temp", "salt", "u", "v")) {
    temp_var <- ncvar_get(temp_nc, ivar)
    
    goa <- sp::SpatialPointsDataFrame(
      coords = cbind(lon = as.numeric(arctic_lon[lon_idx, lat_idx]),
                     lat = as.numeric(arctic_lat[lon_idx, lat_idx])),
      data = data.frame(var = as.numeric(temp_var[lon_idx, 
                                                  lat_idx, 
                                                  depth_layer])))
    
    goa_ras <- raster::raster(x = goa,
                              nrows = diff(range(lon_idx)),
                              ncols = diff(range(lat_idx)))
    goa_ras <- raster::rasterize(x = goa,
                                 y = goa_ras,
                                 field = "var",
                                 na.rm = FALSE)
    
    projection(goa_ras) <- CRS("+proj=longlat +datum=WGS84")
    goa_ras <- projectRaster(goa_ras, crs = crs(AK_land))
    
    assign(value = raster::brick(list(get(paste0(ivar, "_brick")),
                                      goa_ras)), 
           x = paste0(ivar, "_brick"))
    
  }
}

temp_idx <- (length(unique_dates) + 1):length(names(temp_brick))

names(temp_brick)[temp_idx] <- names(v_brick)[temp_idx] <- names(u_brick)[temp_idx] <- names(salt_brick)[temp_idx] <- paste0("X", c(1990, 2012), "_Aug_Avg")

##################################################
####   Plot Covariates
##################################################
# par(mfrow = c(2, 2), mar = c(1, 1, 1, 4))
# for(icovar in c(names(covariate_names_bathy)) ) {
#   plot(static_covariate_brick[[icovar]], axes = F)
#   mtext(text = icovar, side = 3, line = -1.5)
#   plot(AK_land, add = TRUE, col = "tan", border = FALSE)
# }

##################################################
####   Write Rasters
##################################################
for (ifile in  c("static_covariate_brick", 
                 "temp_brick", "salt_brick", "u_brick", "v_brick")) {
  writeRaster(x = get(ifile), 
              filename = paste0("data/covariate_data/", ifile, ".grd"))
}

