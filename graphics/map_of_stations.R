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
# library(RColorBrewer)
# library(plotrix)

##################################################
####   Import AK land
##################################################
AK_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
CA_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/canada_dcw.shp")
RUS_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/russia_dcw.shp")

##################################################
####   Import data
##################################################
otter <- read.csv(paste0("data/fish_data/otter_trawl/",
                         "AK_BTS_Arctic_processed_wide.csv"))
beam <- subset(x = read.csv(paste0("data/fish_data/2017_2019_Beam/",
                                   "ierl_data_processed.csv")),
               common_name == "Arctic cod")

otter_coords <- otter[otter$GEAR_CAT == "otter", 
                      c("YEAR", "MEAN_LONGITUDE", "MEAN_LATITUDE")]
names(otter_coords) <- c("year", "lon", "lat")

beam_coords <- otter[otter$GEAR_CAT == "beam", 
                     c("YEAR", "MEAN_LONGITUDE", "MEAN_LATITUDE")]
names(beam_coords) <- c("year", "lon", "lat")
beam_coords <- rbind(beam_coords, beam[, c("year", "lon", "lat")])  

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
xrange <- extent(otter_sp)[1:2]
yrange <- extent(otter_sp)[3:4]

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
#### Plot Stations
##################################################
{
  png(filename = "graphics/map_of_stations.png", 
      width = 6, height = 3, units = "in", res = 500)
  
  par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
  for(iyear in unique(beam_coords$year)) {
    plot(1, type = "n", asp = 1, axes = F, ann = F,
         xlim = xrange, ylim = yrange)
    points(beam_sp[beam_coords$year == iyear], pch = 16)
    plot(AK_land, add = TRUE, col = "tan", border = FALSE)
    plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
    
    text(x = xrange[1] + diff(xrange) * 1,
         y = yrange[1] + diff(yrange) * 0.25,
         labels = paste0("Beam Trawl\n", iyear))
    box()
  }
  for(iyear in sort(unique(otter_coords$year))) {
    plot(1, type = "n", asp = 1, axes = F, ann = F,
         xlim = xrange, ylim = yrange)
    points(otter_sp[otter_coords$year == iyear], pch = 16)
    plot(AK_land, add = TRUE, col = "tan", border = FALSE)
    plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
    
    text(x = xrange[1] + diff(xrange) * 1,
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


plot(AK_land)
plot(SpPolygon)

