###############################################################################
## Project:      Show location of beam/otter trawl stations 
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(rgdal)
library(sp)
library(raster)
library(graticule)
library(readxl)
# library(RColorBrewer)
# library(plotrix)

##################################################
####   Import AK land
##################################################
AK_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
CA_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/canada_dcw.shp")
RUS_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/russia_dcw.shp")

##################################################
####   Import Chukchi Grid
##################################################
chukchi_grid <- rgdal::readOGR(paste0("data/spatial_data",
                                      "/survey_boundary/CHUKCHI_2012.shp"))

##################################################
####   Import data
##################################################
otter <- read.csv(paste0("data/fish_data/otter_trawl/",
                         "AK_BTS_Arctic_processed_wide.csv"))
beam <- subset(x = read.csv(paste0("data/fish_data/2017_2019_Beam/",
                                   "ierl_data_processed.csv")),
               common_name == "Arctic cod")

beam_norcross <- 
  as.data.frame(read_xlsx(paste0("data/fish_data/Norcross_Beam/",
                                 "2004-09 PSBT fish data Chukchi",
                                 "_for SOAR_18Mar14.xlsx"), 
                          sheet = "static CATCH per haul", 
                          skip = 3))
beam_norcross$lon <- ifelse(beam_norcross$LongitudeStart < 0, 
                            beam_norcross$LongitudeStart, 
                            -1 * beam_norcross$LongitudeStart) 
beam_norcross$lat <- beam_norcross$LatitudeStart
beam_norcross$year <- as.numeric(substr(x = beam_norcross$HaulUniqueCDF, 
                                        start = 1, stop = 4))

otter_coords <- otter[otter$GEAR_CAT == "otter", 
                      c("YEAR", "MEAN_LONGITUDE", "MEAN_LATITUDE")]
names(otter_coords) <- c("year", "lon", "lat")

beam_coords <- otter[otter$GEAR_CAT == "beam", 
                     c("YEAR", "MEAN_LONGITUDE", "MEAN_LATITUDE")]
names(beam_coords) <- c("year", "lon", "lat")
beam_coords <- rbind(beam_coords, 
                     beam[, c("year", "lon", "lat")],
                     beam_norcross[, c("year", "lon", "lat")])  

##################################################
####   Create spatial datapoints
##################################################
beam_sp <- sp::SpatialPoints(coords = beam_coords[, c("lon", "lat")], 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))
beam_sp <- sp::spTransform(x = beam_sp, CRSobj = raster::crs(AK_land))

otter_sp <- sp::SpatialPoints(coords = otter_coords[, c("lon", "lat")], 
                              proj4string = CRS("+proj=longlat +datum=WGS84"))
otter_sp <- sp::spTransform(x = otter_sp, CRSobj = raster::crs(AK_land))

##################################################
#### Map inset snapshot box
##################################################
xrange <- extent(chukchi_grid)[1:2]
yrange <- extent(chukchi_grid)[3:4]

x1 = min(otter$MEAN_LONGITUDE)
x2 = max(otter$MEAN_LONGITUDE)
y1 = min(otter$MEAN_LATITUDE)
y2 = max(otter$MEAN_LATITUDE)

myPolygon = Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))

myPolygons = Polygons(list(myPolygon), ID = "A")
SpPolygon = SpatialPolygons(list(myPolygons), 
                            proj4string = crs("+proj=longlat +datum=WGS84"))
SpPolygon =sp::spTransform(x = SpPolygon, CRSobj = crs(AK_land))

##################################################
#### Plot Beam Trawl Stations
##################################################
{
  png(filename = "graphics/map_of_beamtrawl_stations.png", 
      width = 5, height = 5, units = "in", res = 500)
  
  par(mfrow = c(3, 3), mar = c(0, 0, 0, 0))
  for(iyear in sort(unique(beam_coords$year))) {
    plot(1, type = "n", asp = 1, axes = F, ann = F,
         xlim = xrange, ylim = yrange)
    points(beam_sp[beam_coords$year == iyear], pch = 16)
    plot(AK_land, add = TRUE, col = "tan", border = FALSE)
    plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
    lines(chukchi_grid)
    
    text(x = xrange[1] + diff(xrange) * 0.8,
         y = yrange[1] + diff(yrange) * 0.25,
         labels = paste0("Beam Trawl\n", iyear))
    box()
  }
  
  ##################################################
  #### Map Inset
  ##################################################
  ## Plot legend
  par(mar = c(0, 1, 2, 1))
  plot(AK_land, col = "tan", border = F)
  plot(CA_land, col = "tan", border = F, add = TRUE)
  plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
  # box()
  
  lons <- seq(from = -130, to = -170, by = -10)
  lats <- seq(from = 55, to = 70, by = 5)
  # optionally, specify the extents of the meridians and parallels
  # here we push them out a little on each side
  xl <-  range(lons) + c(0.5, -0.4)
  yl <- range(lats) + c(0.5, -0.4)
  
  grat <- graticule(lons, lats, proj = crs(AK_land))
  grat_labs <- graticule(lons, lats, 
                         xlim = xl, ylim = yl,
                         proj = crs(AK_land), )
  labs <- graticule_labels(lons, lats, 
                           xline = min(xl), yline = min(yl), 
                           proj = crs(AK_land))
  plot(grat, add = TRUE, col = "darkgrey")
  
  text(subset(labs, labs$islon)[c(2:4), ],
       lab = parse(text = labs$lab[labs$islon][c(2:4)]),
       pos = 1, xpd = TRUE, cex = 0.7)
  text(subset(labs, !labs$islon),
       lab = parse(text = labs$lab[!labs$islon]),
       pos = 2, xpd = TRUE, cex = 0.7)
  
  plot(SpPolygon, add = TRUE, lwd = 2)
  
  dev.off()
}

##################################################
#### Plot Otter Trawl Stations
##################################################
{
  png(filename = "graphics/map_of_ottertrawl_stations.png", 
      width = 5, height = 2, units = "in", res = 500)
  
  par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))
  for(iyear in sort(unique(otter_coords$year))) {
    plot(1, type = "n", asp = 1, axes = F, ann = F,
         xlim = xrange, ylim = yrange)
    points(otter_sp[otter_coords$year == iyear], pch = 16)
    plot(AK_land, add = TRUE, col = "tan", border = FALSE)
    plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
    lines(chukchi_grid)
    
    text(x = xrange[1] + diff(xrange) * 0.8,
         y = yrange[1] + diff(yrange) * 0.25,
         labels = paste0("Otter Trawl\n", iyear))
    box()
  }
  
  ##################################################
  #### Map Inset
  ##################################################
  ## Plot legend
  par(mar = c(0, 1, 2, 1))
  plot(AK_land, col = "tan", border = F)
  plot(CA_land, col = "tan", border = F, add = TRUE)
  plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
  # box()
  
  lons <- seq(from = -130, to = -170, by = -10)
  lats <- seq(from = 55, to = 70, by = 5)
  # optionally, specify the extents of the meridians and parallels
  # here we push them out a little on each side
  xl <-  range(lons) + c(0.5, -0.4)
  yl <- range(lats) + c(0.5, -0.4)
  
  grat <- graticule(lons, lats, proj = crs(AK_land))
  grat_labs <- graticule(lons, lats, 
                         xlim = xl, ylim = yl,
                         proj = crs(AK_land), )
  labs <- graticule_labels(lons, lats, 
                           xline = min(xl), yline = min(yl), 
                           proj = crs(AK_land))
  plot(grat, add = TRUE, col = "darkgrey")
  
  text(subset(labs, labs$islon)[c(2:4), ],
       lab = parse(text = labs$lab[labs$islon][c(2:4)]),
       pos = 1, xpd = TRUE, cex = 0.7)
  text(subset(labs, !labs$islon),
       lab = parse(text = labs$lab[!labs$islon]),
       pos = 2, xpd = TRUE, cex = 0.7)
  
  plot(SpPolygon, add = TRUE, lwd = 2)
  
  dev.off()
}



