###############################################################################
## Project:         Fiddling with Arctic ROMS output
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################

##################################################
#### Import Packages
##################################################
library(ncdf4)

##################################################
#### Import the arctic grid in the PAROMS folder as well as an example nc
#### file from January 1st 1990.
##################################################
LB_dir <- c("G:/Barnett/PAROMS/Arctic4/")

arctic_grid <- ncdf4::nc_open(paste0(LB_dir, "grid_Arctic_4.nc"))
example_nc <- ncdf4::nc_open(paste0(LB_dir, 
                                    "averages/arctic4_avg_1990-01-01T.nc"))

##################################################
#### truncated output from print(example_nc) for the temperature and salt variables
##################################################
# float temp[xi_rho,eta_rho,s_rho,ocean_time]   (Chunking: [173,273,13,1])  (Compression: shuffle,level 1)
# long_name: time-averaged potential temperature
# units: Celsius
# time: ocean_time
# grid: grid
# location: face
# coordinates: lon_rho lat_rho s_rho ocean_time
# field: temperature, scalar, series
# _FillValue: 9.99999993381581e+36

# float salt[xi_rho,eta_rho,s_rho,ocean_time]   (Chunking: [173,273,13,1])  (Compression: shuffle,level 1)
# long_name: time-averaged salinity
# time: ocean_time
# grid: grid
# location: face
# coordinates: lon_rho lat_rho s_rho ocean_time
# field: salinity, scalar, series
# _FillValue: 9.99999993381581e+36

##################################################
#### temperature and salt have dims lon_rho, lat_rho, s_rho, and ocean_time
#### the ocean_time variable is of length 1 (daily average). The rho refers to 
#### its position on the arakawa C-grid 
#### (https://www.myroms.org/wiki/Fractional_Coordinate_System_(%CE%BE_-_%CE%B7_space))
#### s_rho I imagine refers to depth???
##################################################

##################################################
#### I can find lat and lon values from the arctic grid
##################################################
lon <- matrix(ncdf4::ncvar_get(arctic_grid, "lon_rho"))
lat <- matrix(ncdf4::ncvar_get(arctic_grid, "lat_rho"))

locs_df <- data.frame(lon = unlist(lon),
                      lat = unlist(lat))

##################################################
#### extract temperature data
##################################################
temperature <- ncdf4::ncvar_get(nc = example_nc, varid = "temp")
temperature_bottom <- temperature[, , 1] 

temperature_bottom_df <- data.frame(temp = unlist(matrix(temperature_bottom)))

locs_df$temp <- temperature_bottom_df$temp + 1 + abs(min(temperature_bottom_df$temp, na.rm= TRUE))


##################################################
#### 
##################################################



##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

##################################################
#### 
##################################################

