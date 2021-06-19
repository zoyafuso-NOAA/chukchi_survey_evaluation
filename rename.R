###############################################################################
## Project:       VAST Modelling, Chukchi Sea
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   1990 and 2012 Chukchi data
###############################################################################

library(rgdal)
library(sp)

##################################################
####   Import Data
##################################################
ak_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
extrapolation_grid <- read.csv("data/spatial_data/BS_Chukchi_extrapolation_grids/ChukchiThorsonGrid.csv")

##################################################
####   Synthesize CVs
##################################################
cvs <- data.frame()

for (igear in c("beam_trawl", "otter_trawl")) {
  spp_list <- gsub(x = list.dirs(path = paste0("results/", igear), 
                                 recursive = F),
                   pattern = paste0("results/", igear, "/"),
                   replacement = "")
  
  years_to_include <- list(beam_trawl = c(2012, 2017, 2019),
                           otter_trawl = c(1990, 2012))[[igear]]
  
  for (ispp in spp_list) {
    ## Extract CVs
    temp_cvs <- read.csv(paste0("results/", igear, "/", 
                                ispp, "/Table_for_SS3.csv"))
    temp_cvs <- subset(x = temp_cvs, 
                       subset = Year %in% years_to_include)$SD_log
    cvs <- rbind(cvs,
                 data.frame(species = ispp,
                            gear = igear,
                            year = years_to_include,
                            cv = round(temp_cvs, 2)))
  }
}

cv_table <- list(
  otter = tidyr::spread(data = subset(cvs, gear == "otter_trawl"), 
                        key = year, value = cv),
  beam = tidyr::spread(data = subset(cvs, gear == "beam_trawl"), 
                       key = year, value = cv))
save(cv_table, file = "presentations/results_06_19_2021/cv_table.RDS")

##################################################
####   Number of stations
##################################################
BT_data_2017_2019 <- 
  read.csv("data/fish_data/2017_2019_Beam/ierl_data_processed.csv")
BT_data_2017_2019 <- subset(x = BT_data_2017_2019,
                            subset = common_name == "Arctic cod")
table(BT_data_2017_2019$year)

BTS_data <- read.csv("data/fish_data/otter_trawl/AK_BTS_Arctic_processed.csv")
BTS_data <- subset(x = BTS_data,
                   subset = common_name == "Arctic cod")
with(subset(x = BTS_data, subset = common_name == "Arctic cod"),
     table(year, gear))

##################################################
####   map of Stations
##################################################
all_pts <- rbind(BT_data_2017_2019[, c("year", "lat", "lon", "gear")],
                 BTS_data[, c("year", "lat", "lon", "gear")])

all_pts_utm <- sp::SpatialPointsDataFrame(
  coords = all_pts[, c("lon", "lat")],
  proj4string = CRS("+proj=longlat +datum=WGS84"), 
  data = all_pts)
all_pts_utm <- sp::spTransform(x = all_pts_utm,
                               CRSobj = CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 + x_0=0 +y_0=0 +datum=NAD83 +units=m"))

{
  png(filename = c("graphics/map_of_stations.png"), 
      height = 2.5, width = 4, units = "in", res = 500)
  par(mar = c(0, 0, 3, 0), mfrow = c(2, 3), oma = c(0, 0, 0, 0))
  for (igear in c("beam", "otter")) {
    years_to_include <- sort(unique(all_pts$year[all_pts$gear == igear]))
    for (iyear in years_to_include) {
      plot(subset(all_pts_utm, year == 2012 & gear == igear), 
           xlim = all_pts_utm@bbox["lon", ], 
           ylim = all_pts_utm@bbox["lat", ],
           main = paste(igear, "trawl,", iyear),
           cex = 0.5)
      plot(ak_land, add = TRUE, col = "tan", border = F)
      box()
    }
  }
  dev.off()
}

##################################################
####   Plot Densities
##################################################

##################################################
####   Mapping function
##################################################
make_a_raster <- function(extrap_grid = extrapolation_grid,
                          lat_lon_names = c("Lon", "Lat"),
                          plot_what = fit_sub$D_gct[, 1, 1],
                          raster_res = 0.11) {
  
  plot_layer <- sp::SpatialPointsDataFrame(
    coords = extrapolation_grid[, lat_lon_names],
    data = data.frame(plot_this = plot_what) )
  plot_ras <- raster::raster(x = plot_layer, 
                             resolution = raster_res)
  plot_ras <- raster::rasterize(x = plot_layer, 
                                y = plot_ras, 
                                field = "plot_this")
  
  return(plot_ras)
}

spp_list <- sort(unique(gsub(x = list.dirs(path = paste0("results/beam_trawl"), 
                                           recursive = F),
                             pattern = paste0("results/beam_trawl/"),
                             replacement = ""),
                        
                        gsub(x = list.dirs(path = paste0("results/otter_trawl"), 
                                           recursive = F),
                             pattern = paste0("results/otter_trawl/"),
                             replacement = "")
))


for (ispp in spp_list) {
  
  png(filename = paste0(getwd(), "/presentations/results_06_19_2021/", 
                        ispp, "_density.png"),
      width = 5, height = 4, units = "in", res = 500)
  par(mfrow = c(2, 3), mar = c(3, 3, 1, 3), oma = c(0, 0, 3, 1))
  for (igear in c("beam", "otter")) {
    years_to_include <- list(beam = c(1, 6, 8),
                             otter = c(1, 23))[[igear]]
    years_labels <- list(beam = c(2012:2019),
                         otter = c(1990:2012))[[igear]]
    
    # n_spp <- length(spp_list)
    n_year <- length(years_to_include)
    
    file_name <- paste0("results/", igear, "_trawl/", ispp, "/fit.RData")
    if(file.exists(file_name)) {
      load(file = file_name)
      
      for (iyear in years_to_include) {
        plot(make_a_raster(fit_sub$D_gct[, , iyear]))
        mtext(side = 3, text = paste(igear, "trawl,", years_labels[iyear]))
      }
      
    }
    
  }
  
  mtext(side = 3, text = ispp, outer = T, line = 1, font = 2)
  dev.off()
}

