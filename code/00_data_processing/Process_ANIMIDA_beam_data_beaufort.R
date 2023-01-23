###############################################################################
## Project:       Synthesize Alaska Beam Trawl Beaufort Data from 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries
##################################################
library(reshape)
library(rgdal)
library(raster)

##################################################
####  Import raw data
##################################################
data <- (read.csv("data/fish_data/2014 Beaufort beam/0162530/1.1/data/0-data/ANIMIDA-III-Sample-Results/CSV/DataValues.csv"))
taxa <- (read.csv("data/fish_data/2014 Beaufort beam/0162530/1.1/data/0-data/ANIMIDA-III-Sample-Results/CSV/Taxonomy.csv"))
depth_ras <- raster::raster(
  x = "data/covariate_data/Arctic Bathymetry 1 km_NEW/arctic1km_new/")

##################################################
####  Species list (Arctic cod is the only one on the list present)
##################################################
spp_list <- c("Arctic cod" = "Boreogadus saida", 
              "saffron cod" = "Eleginus gracilis",
              "walleye pollock" =  "Gadus chalcogrammus",
              "Pacific cod" = "Gadus macrocephalus",
              "Bering flounder" = "Hippoglossoides robustus",
              "yellowfin sole" = "Limanda aspera",
              "Alaska plaice" = "Pleuronectes quadrituberculatus")
taxa_ids <- taxa$TaxaID[taxa$Species %in% spp_list]

##################################################
####  Subset data to only abundance, bottom trawl data for the species in 
####     spp_list. Select relevant fields. Convert individuals / 1000 m^2 to 
####     individuals / km^2 by multiplying by 1000
##################################################
cpue <- subset(x = data, 
               subset = SampleType == "Trawl" & 
                 SampleMedium == "Benthic zone" & 
                 VariableName == "Abundance" & 
                 TaxaID %in% taxa_ids,
               select = c(DataValue, Latitude, Longitude, LocalDateTime) )
cpue$year <- as.numeric(substr(x = cpue$LocalDateTime, start = 1, stop = 4))
cpue$spp_name <- "Arctic cod"
cpue$DataValue <- cpue$DataValue* 1000

##################################################
####  Rename fields
##################################################
names(cpue) <- c("ind_km2", "lat", "lon", "date", "year", "spp_name")

##################################################
####  Extract depths from Pirtle raster
##################################################
cpue_pts <- 
  sp::SpatialPoints(coords = cpue[, c("lon", "lat")], 
                    proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
cpue_pts_aea <- sp::spTransform(x = cpue_pts,
                                CRSobj = raster::crs(depth_ras) )

cpue$depth <- -1 * raster::extract(x = depth_ras, y = cpue_pts_aea)
cpue$gear <- "beam"

##################################################
####  Create interpolation grid: create buffer around sp_poly
##################################################
coord_input <- cpue[, c("lon", "lat")]
names(coord_input) <- c("Lon", "Lat")

beaufort_input <- FishStatsUtils::make_extrapolation_info(
  Region = "other",
  grid_dim_km = c(3.704, 3.704), #2 nm
  observations_LL = coord_input)

beaufort_grid <- beaufort_input$Data_Extrap

grid_pts <- 
  sp::SpatialPoints(coords = beaufort_grid[, c("Lon", "Lat")],
                    proj4string = sp::CRS("+proj=longlat +datum=WGS84"))

grid_pts_aea <- sp::spTransform(x = grid_pts, 
                                CRSobj = raster::crs(depth_ras))

beaufort_grid$depth <- raster::extract(x = depth_ras, y = grid_pts_aea) * - 1
beaufort_grid <- subset(x = beaufort_grid, 
                        subset = depth > floor(min(cpue$depth)) &
                          depth < ceiling(max(cpue$depth)))

##################################################
####  Save
##################################################
write.csv(x = cpue, 
          file = paste0("data/fish_data/2014 Beaufort beam/",
                        "Beaufort_beam_processed_long.csv"),
          row.names = FALSE)

write.csv(x = beaufort_grid, 
          file = paste0("data/spatial_data/Beaufort_extrapolation_grids/",
                        "Beaufort_Grid.csv"),
          row.names = FALSE)
