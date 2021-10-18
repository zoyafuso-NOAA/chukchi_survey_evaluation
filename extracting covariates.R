rm(list = ls())

library(ncdf4); library(raster); library(rgdal)
akland <- readOGR("G:/Oyafuso/Arctic/spatial_data/land_shapefiles/AKland.shp")
covariate_dir <- "G:/Barnett/PAROMS/Arctic4/averages/"

otter_df <- read.csv(file = "data/fish_data/otter_trawl/AK_BTS_Arctic_processed_wide.csv")

arctic_grid <- ncdf4::nc_open(paste0(dirname(covariate_dir), "/grid_Arctic_4.nc"))

arctic_lon <- ncvar_get(arctic_grid, "lon_rho")
arctic_lon <- ifelse(arctic_lon > 180, -360 + arctic_lon, arctic_lon)
summary(as.numeric(unlist(arctic_lon)))

arctic_lat <- ncvar_get(arctic_grid, "lat_rho")

plot(unlist(arctic_lon), unlist(arctic_lat))

unique_dates <- unique(otter_df$DATE)

for (idate in unique_dates[1]){
  iyear <- substr(x = idate, start = 1, stop = 4)
  imonth <- substr(x = idate, start = 6, stop = 7)
  iday <- substr(x = idate, start = 9, stop = 10)
  
  ifile <- paste0(covariate_dir, 
                  "arctic4_avg_", iyear, "-", imonth, "-", iday, "T.nc")
  
  temp_nc <- ncdf4::nc_open(filename = ifile)
  temperature <- ncvar_get(temp_nc, "temp")
  salt <- ncvar_get(temp_nc, "salt")
  x_curr <- ncvar_get(temp_nc, "u")
  y_curr <- ncvar_get(temp_nc, "v")
}

lon_idx <- 310:600; lat_idx <- 51:550; depth_layer <- 1

goa <- sp::SpatialPointsDataFrame(
  coords = cbind(lon = as.numeric(arctic_lon[lon_idx, lat_idx]),
                 lat = as.numeric(arctic_lat[lon_idx, lat_idx])),
  data = data.frame(temperature = as.numeric(temperature[lon_idx, lat_idx, depth_layer]),
                    salt = as.numeric(salt[lon_idx, lat_idx, depth_layer]),
                    x_curr = as.numeric(x_curr[lon_idx, lat_idx, depth_layer]),
                    y_curr = as.numeric(y_curr[lon_idx, lat_idx, depth_layer])) )

goa_ras <- raster::raster(x = goa,
                          nrows = diff(range(lon_idx)),
                          ncols = diff(range(lat_idx)))
goa_ras <- raster::rasterize(x = goa,
                             y = goa_ras,
                             field = "y_curr",
                             na.rm = FALSE)
projection(goa_ras) <- CRS("+proj=longlat +datum=WGS84")

# plot(raster::projectRaster(goa_ras, crs = "+proj=utm +zone=5N +units=km"))
# plot(goa_ras)

par(mar = c(1, 1, 3, 3))
plot(projectRaster(goa_ras, crs = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "), axes = F,
     main = "y-current")
plot(akland, add = TRUE, col = "brown")
