##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       NBS bottom trawl stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Lewis Barnett (lewis.barnett@noaa.gov)
## Description:   Constants used in survey optimization and simulations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(terra)
library(FishStatsUtils)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Shapefiles   ----
##   Save the crs for the lat/lon projection for future use
##   Save the crs used for AKland to have a consistent Albers projection
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AKland <- terra::vect("data/spatial_data/land_shapefiles/AKland.shp")
aea_crs <- terra::crs(AKland)
latlon_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Dist to shore ----
##   Calculate minimum distance to land for each grid point
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
grid_pts <- 
  terra::vect(x = FishStatsUtils::northern_bering_sea_grid[, c("Lon", "Lat")], 
              type = "points",
              crs = latlon_crs)
grid_pts <- terra::project(x = grid_pts,  aea_crs)

## Calculate the minimium distance between each point 
## and the many polygons in AKland
dist_to_shore <- apply(X = terra::distance(x = grid_pts, y = AKland, 
                                           pairwise = FALSE), 
                       MARGIN = 1, FUN = min)

## variables for plotting 
xrange <- terra::ext(grid_pts)[1:2]
yrange <- terra::ext(grid_pts)[3:4]

n_cells <- nrow(FishStatsUtils::northern_bering_sea_grid)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species Included for each gear ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_list <- read.csv(file = "results/good_species.csv")
n_spp <- nrow(spp_list)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Years sampled ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
years <- c(1985, 1988, 1991, 2010, 2017:2019, 2021:2022)
n_years <- length(years)
year_idx <- which(1985:2022 %in% years)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Number of simulations ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_iters <- 1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   SamplingStrata Data Input ----
##   X1 and X2 are stratum variables latitude (X1) and distance from shore (X2)
##   Y1, Y2, ..., Yn_spp are density predictions from VAST and 
##   Y1_SUM_SQ, ..., is the square of the density prediction
##
##   true_index is the true index calculated from VAST
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
frame_df<- data.frame(domainvalue = 1,
                      id = 1:nrow(FishStatsUtils::northern_bering_sea_grid),
                      X1 = FishStatsUtils::northern_bering_sea_grid[, "Lat"],
                      X2 = dist_to_shore)

true_index <- matrix(nrow = length(x = year_idx), ncol = n_spp,
                     dimnames = list(paste0("Year_", years), 
                                     paste0("Y", 1:n_spp)))

frame_df$WEIGHT <- length(x = years)
idir <- paste0("nbs_otter/vast_fits/")

## Extract predicted densities and indces for each taxon
for (ispp in 1:n_spp) { ## Loop over taxon -- start
  load(file = paste0("results/", idir, spp_list$taxon[ispp], "/fit.RData"))
  frame_df[, paste0("Y", ispp)] <- 
    rowSums(x = as.matrix(x = fit$Report$D_gct[, 1, year_idx]))
  frame_df[, paste0("Y", ispp, "_SQ_SUM")] <- 
    rowSums(x = as.matrix(x = fit$Report$D_gct[, 1, year_idx]^2))
  
  true_index[, ispp] <- fit$Report$Index_ctl[1, year_idx, 1]
} ## Loop over taxon -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
save(list = c("aea_crs", "latlon_crs", "xrange", "yrange", 
              "n_cells", "n_iters", "spp_list", "n_spp", 
              "n_years", "frame_df", "true_index"), 
     file = "data/survey_opt_data/optimization_data_nbs.RData")

terra::writeVector(x = grid_pts, file = "data/survey_opt_data/grid_pts_nbs.shp",
                   overwrite = TRUE)

