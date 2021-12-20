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
n_iters <- 500

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
spp_list_otter <- c("Alaska plaice", "Arctic cod", "Bering flounder", 
                    "saffron cod", "snow crab", "yellowfin sole")
n_spp_otter <- length(spp_list_otter)
names(spp_list_otter) <- paste0("Y", 1:n_spp_otter)

spp_list_beam <- c("Arctic cod", "Bering flounder", 
                   "saffron cod", "snow crab", "yellowfin sole")
n_spp_beam <- length(spp_list_beam)
names(spp_list_beam) <- paste0("Y", 1:n_spp_beam)

##################################################
#### Constants related to gear type
##################################################
years_otter <- c(1990, 2012)
year_idx_otter <- c(1, 23)

years_beam <- c(2012, 2017, 2019)
year_idx_beam <- c(1, 6, 8)

##################################################
#### Species optimization data input
#### X1 and X2 are stratum variables latitude (X1) and distance from shore (X2)
#### Y1, Y2, ..., Yn_spp are density predictions from VAST and 
#### Y1_SUM_SQ, ..., is the square of the density prediction
####
#### true_index is the true index from VAST
##################################################
frame_df_otter <- frame_df_beam <- data.frame(domainvalue = 1,
                                              id = 1:nrow(chukchi_grid),
                                              X1 = chukchi_grid$Lat,
                                              X2 = dist_to_shore)

for (igear in c("otter", "beam")) {
  
  spp_list <- get(paste0("spp_list_", igear))
  n_spp <- get(paste0("n_spp_", igear))
  
  years <- get(paste0("years_", igear))
  year_idx <- get(paste0("year_idx_", igear))

  true_index <- matrix(nrow = length(year_idx), ncol = n_spp,
                       dimnames = list(paste0("Year_", years), 
                                       paste0("Y", 1:n_spp)))
  
  frame_df <- get(paste0("frame_df_", igear))
  frame_df$WEIGHT <- length(years)
  idir <- ifelse(test = igear == "otter", 
                 yes = "otter_trawl/vast_fits/", 
                 no = "beam_trawl_2012_2019/vast_fits/")

  
  for (ispp in 1:n_spp) {
    load(paste0("results/", idir, spp_list[ispp], "/fit.RData"))
    frame_df[, paste0("Y", ispp)] <- 
      rowSums(as.matrix(fit$Report$D_gct[, 1, year_idx]))
    frame_df[, paste0("Y", ispp, "_SQ_SUM")] <- 
      rowSums(as.matrix(fit$Report$D_gct[, 1, year_idx]^2))
    
    true_index[, ispp] <- fit$Report$Index_ctl[1, year_idx, 1]
  }
  
  assign(x = paste0("frame_df_", igear), value = frame_df)
  assign(x = paste0("true_index_", igear), value = true_index)  
  
  rm(igear, spp_list, n_spp, years, year_idx, true_index, frame_df, idir, ispp)
}

##################################################
#### Save
##################################################
save(list = c("aea_crs", "latlon_crs", "xrange", "yrange", "grid_pts", 
              "n_cells", "n_iters", 
              as.vector(sapply(X = c( "spp_list", "n_spp", 
                                     "frame_df", "true_index"), 
                               FUN = function(x) 
                                 paste0(x, c("_beam", "_otter"))))),
     file = "data/survey_opt_data/optimization_data.RData")

