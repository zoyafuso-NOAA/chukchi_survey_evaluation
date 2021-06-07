###############################################################################
## Project:      
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  
###############################################################################
rm(list = ls())

library(rgdal)
library(sp)
library(raster)

data_dir <- "G:/Oyafuso/Arctic/spatial_data/"

nbs_grid <- read.csv( file = paste0(data_dir, 
                                    "/BS_Chukchi_extrapolation_grids",
                                    "/NBSThorsonGrid.csv") )
chukchi_grid <- read.csv( file = paste0(data_dir, 
                                        "/BS_Chukchi_extrapolation_grids",
                                        "/ChukchiThorsonGrid.csv") )
ebs_grid <- read.csv( file = paste0(data_dir, 
                                    "/BS_Chukchi_extrapolation_grids",
                                    "/EBSThorsonGrid.csv") )

bs_chukchi_grid <- rbind(cbind(Region = "EBS",
                               ebs_grid),
                         cbind(Region = "NBS",
                               nbs_grid),
                         cbind(Region = "CHUKCHI",
                               chukchi_grid))

AK_land <- rgdal::readOGR(paste0(data_dir, "/land_shapefiles/AKland.shp")) 


grid_pts_latlon <- sp::SpatialPoints(
  coords = bs_chukchi_grid[, c("ShapeLon", "ShapeLat")],
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

grid_pts_aea <- sp::spTransform(x = grid_pts_latlon,
                                CRSobj = CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 + x_0=0 +y_0=0 +datum=NAD83 +units=m") )

ebs_boundary <- raster::shapefile(paste0(data_dir, 
                                         "/survey_boundary/EBS_NBS_2019.shp"))
chukchi_boundary <- raster::shapefile(paste0(data_dir, 
                                             "/survey_boundary/CHUKCHI_2012.shp"))


par(mar = c(0, 0, 0, 0))
plot(grid_pts_aea@coords,
     pch = 16,
     cex = 0.05,
     asp = 1,
     axes = F,
     col = "white")
plot(ebs_boundary, add = T, col = "black")
plot(chukchi_boundary, add = T, col = "cyan", border = F)
plot(AK_land, add = T, col = "tan", border = "tan")

text(x = c(mean(ebs_boundary@bbox["x", ]) - 150000,
           mean(chukchi_boundary@bbox["x", ]) + 50000,
           mean(chukchi_boundary@bbox["x", ]) + 600000),
     y = c(mean(ebs_boundary@bbox["y", ]),
           mean(chukchi_boundary@bbox["y", ]) + 200000,
           mean(chukchi_boundary@bbox["y", ]) + 200000),
     labels = c("NBS and \n EBS", "Chukchi", "Beaufort"),
     col = c("white", "black", "black"), 
     cex = 2)

