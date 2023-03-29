##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Chukchi bottom trawl stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
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
  terra::vect(x = FishStatsUtils::chukchi_sea_grid[, c("Lon", "Lat")], 
              type = "points",
              crs = latlon_crs)
grid_pts <- terra::project(x = grid_pts,  aea_crs)

## Calculate the minimium distance between each point 
## and the many polygons in AKland
dist_to_shore <- apply(X = terra::distance(x = grid_pts, y = AKland, 
                                           pairwise = FALSE, unit = "km"), 
                       MARGIN = 1, FUN = min)

## variables for plotting 
xrange <- terra::ext(grid_pts)[1:2]
yrange <- terra::ext(grid_pts)[3:4]

n_cells <- nrow(FishStatsUtils::chukchi_sea_grid)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species Included for each gear ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_list <- read.csv(file = "results/good_species.csv")

spp_list_otter <- spp_list$taxon[spp_list$otter]
n_spp_otter <- length(spp_list_otter)
names(spp_list_otter) <- paste0("Y", 1:n_spp_otter)

spp_list_beam <- spp_list$taxon[spp_list$beam]
n_spp_beam <- length(spp_list_beam)
names(spp_list_beam) <- paste0("Y", 1:n_spp_beam)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Years sampled for each gear type ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
years_otter <- c(1990, 2012)
year_idx_otter <- c(1, 23)

years_beam <- c(2012, 2017, 2019)
year_idx_beam <- c(1, 6, 8)

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
##   true_index is the true index claculated from VAST
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
frame_df_otter <- frame_df_beam <- data.frame(domainvalue = 1,
                                              id = 1:nrow(chukchi_sea_grid),
                                              X1 = chukchi_sea_grid[, "Lat"],
                                              X2 = dist_to_shore)

for (igear in c("otter", "beam")) { ## Loop over gear -- start
  
  ## Grab species list and years based on gear igear
  spp_list <- get(x = paste0("spp_list_", igear))
  n_spp <- get(x = paste0("n_spp_", igear))
  
  years <- get(x = paste0("years_", igear))
  year_idx <- get(x = paste0("year_idx_", igear))
  
  true_index <- matrix(nrow = length(x = year_idx), ncol = n_spp,
                       dimnames = list(paste0("Year_", years), 
                                       paste0("Y", 1:n_spp)))
  
  frame_df <- get(x = paste0("frame_df_", igear))
  frame_df$WEIGHT <- length(x = years)
  idir <- paste0("chukchi_", igear, "/vast_fits/")
  
  ## Extract predicted densities and indces for each taxon
  for (ispp in 1:n_spp) { ## Loop over taxon -- start
    load(file = paste0("results/", idir, spp_list[ispp], "/fit.RData"))
    frame_df[, paste0("Y", ispp)] <- 
      rowSums(x = as.matrix(x = fit$Report$D_gct[, 1, year_idx]))
    frame_df[, paste0("Y", ispp, "_SQ_SUM")] <- 
      rowSums(x = as.matrix(x = fit$Report$D_gct[, 1, year_idx]^2))
    
    true_index[, ispp] <- fit$Report$Index_ctl[1, year_idx, 1]
  } ## Loop over taxon -- end
  
  assign(x = paste0("frame_df_", igear), value = frame_df)
  assign(x = paste0("true_index_", igear), value = true_index)  
  
  rm(igear, spp_list, n_spp, years, year_idx, true_index, frame_df, idir, ispp)
}  ## Loop over gear -- end


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
save(list = c("aea_crs", "latlon_crs", "xrange", "yrange", 
              "n_cells", "n_iters", 
              as.vector(sapply(X = c( "spp_list", "n_spp", 
                                      "frame_df", "true_index"), 
                               FUN = function(x) 
                                 paste0(x, c("_beam", "_otter"))))),
     file = "data/survey_opt_data/optimization_data.RData")


terra::writeVector(x = grid_pts, file = "data/survey_opt_data/grid_pts.shp",
                   overwrite = TRUE)

