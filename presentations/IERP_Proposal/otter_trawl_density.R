###############################################################################
## Project:      VAST density output for the Arctic IERP proposal
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## POC:          Lewis Barnett (lewis.barnett@noaa.gov)
## Description:  Predicted densities by species and year
###############################################################################

rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(rgdal)
library(sp)
library(raster)
library(graticule)
library(RColorBrewer)
library(plotrix)

##################################################
####   Import AK land
##################################################
AK_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
CA_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/canada_dcw.shp")
RUS_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/russia_dcw.shp")

##################################################
####   Import Extrapolation Grid
##################################################
extrapolation_grid <- read.csv(paste0("data/spatial_data/",
                                      "BS_Chukchi_extrapolation_grids/",
                                      "ChukchiThorsonGrid.csv"))
grid_pts_latlon <- sp::SpatialPoints(
  coords = extrapolation_grid[, c("Lon", "Lat")],
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

grid_pts_aea <- sp::spTransform(x = grid_pts_latlon,
                                CRSobj = CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 + x_0=0 +y_0=0 +datum=NAD83 +units=m") )

##################################################
#### Map inset snapshot box
##################################################
x1 = min(extrapolation_grid$Lon)
x2 = max(extrapolation_grid$Lon)
y1 = min(extrapolation_grid$Lat)
y2 = max(extrapolation_grid$Lat)

myPolygon = Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))

myPolygons = Polygons(list(myPolygon), ID = "A")
SpPolygon = SpatialPolygons(list(myPolygons), 
                            proj4string = crs("+proj=longlat +datum=WGS84"))
SpPolygon = sp::spTransform(x = SpPolygon, CRSobj = crs(AK_land))

##################################################
####   Species to plot
##################################################
spp_list <- c("Alaska plaice", "Arctic cod", "Bering flounder", #"saffron cod",
              "snow crab", "walleye pollock" , "yellowfin sole")
years <- c(1990, 2012)

##################################################
####   Plot
##################################################
{
  ## Start device
  png("presentations/IERP_Proposal/otter_trawl_density.png",
      width = 70, height = 200, units = "mm", res = 500)
  
  layout(mat = matrix(c(1:18, 19, 20, 20), ncol = 3, byrow = TRUE),
         widths = c(1, 1, 0.75))
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 5, 0))  
  
  for(ispp in spp_list) {
    load(paste0("results/otter_trawl/", ispp, "/fit.RData"))
    temp_density <- fit$Report$D_gct[, 1, ]
    
    ymax <- max(temp_density)
    
    for (iyear in 1:2) {
      dens_spdf <- sp::SpatialPointsDataFrame(
        coords = grid_pts_aea@coords[, c("Lon", "Lat")],
        data = data.frame(Str_no = temp_density[, iyear]),
        proj4string = crs(grid_pts_aea) )
      
      dens_ras <- raster::raster(x = dens_spdf, resolution = 10000)
      dens_ras <- raster::rasterize(x = dens_spdf,
                                    y = dens_ras,
                                    field = "Str_no")
      
      image(dens_ras, asp = 1, axes = FALSE, zlim = c(0, ymax), col = blues9)

      if(ispp == spp_list[1]) mtext(side = 3,
                                    text = c(1990, 2012)[iyear],
                                    font = 2)
      
      plot(AK_land, add = TRUE, col = "tan", border = FALSE)
      plot(RUS_land, add = TRUE, col = "tan", border = FALSE)
      box()
      
      text(x = extent(dens_ras)[1] + 0.85 * diff(extent(dens_ras)[1:2]),
           y = extent(dens_ras)[3] + 0.30 * diff(extent(dens_ras)[3:4]),
           labels = gsub(x = ispp, replacement = "\n", pattern = " "))
      
      obs_cpue <- sp::SpatialPointsDataFrame(
        coords = fit$data_frame[c("Lon_i", "Lat_i")],
        data = data.frame(year = fit$data_frame$t_i,
                          cex_val = log10((fit$data_frame$b_i / fit$data_frame$a_i) + 1)),
        proj4string = CRS("+proj=longlat +datum=WGS84") )
      
      points(sp::spTransform(subset(x = obs_cpue,
                                    subset = cex_val == 0 &
                                      year == c(1990, 1991)[iyear]),
                             CRSobj = CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 + x_0=0 +y_0=0 +datum=NAD83 +units=m")),
             cex = 0.5,
             pch = "+",
             col = "black")
      
      points(sp::spTransform(subset(x = obs_cpue,
                                    subset = year == c(1990, 1991)[iyear]),
                             CRSobj = CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 + x_0=0 +y_0=0 +datum=NAD83 +units=m")),
             cex = subset(obs_cpue@data, year == c(1990, 1991)[iyear])$cex_val/2,
             pch = 1,
             lwd = 0.5,
             col = "black")
    }
    
    plot(1, type = "n", axes = F, ann = F,
         xlim = c(0, 3), ylim = c(0, 10))
    
    ## Legend
    r.range <- pretty(c(0, ymax), n = 3)
    plotrix::color.legend(
      xl = 0.2,
      xr = 0.7,
      yb = 1,
      yt = 9,
      legend = r.range,
      rect.col = blues9,
      gradient = "y", 
      align = "rb",
      cex = 0.6,
      xpd = NA)
    
  }
  
  ## Plot points legend
  plot(1, axes = F, ann = F, type = "n",
       ylim = c(0, 1), xlim = c(0, 6.5))
  points(y = 0.8, x = 0.25, cex = 1, pch = "+")
  points(y = rep(0.8, 4), x = seq(from = 1.75, by = 1.5, length = 4),
         cex = c((1:4)/2), pch = 1)
  text(y = rep(0.7, 5), x = seq(from = 0.25, by = 1.5, length = 5),
       labels = c(0, 10^(0:3)), pos = 1, xpd = NA, cex = 0.9)
  mtext(side = 1, line = -2.5, text = "Observed CPUE\n(kg/km2)", cex = 0.75)

  ## Plot legend
  par(mar = c(0, 1, 1, 3))
  plot(AK_land, col = "tan", border = F, xlim = c(extent(AK_land)[1] * 0.8,
                                                  extent(AK_land)[2]))
  plot(CA_land, col = "tan", border = F, add = TRUE)
  plot(RUS_land, add = TRUE, col = "tan", border = FALSE)

  lons <- seq(from = -130, to = -170, by = -10)
  lats <- seq(from = 55, to = 70, by = 5)
  # optionally, specify the extents of the meridians and parallels
  # here we push them out a little on each side
  xl <-  range(lons) + c(0.5, -0.4)
  yl <- range(lats) + c(0.5, -0.4)

  grat <- graticule(lons, lats, proj = crs(grid_pts_aea))
  grat_labs <- graticule(lons, lats, proj = crs(grid_pts_aea), xlim = xl, ylim = yl)
  labs <- graticule_labels(lons, lats, xline = min(xl), yline = min(yl), proj = crs(grid_pts_aea))
  plot(grat, add = TRUE, col = "grey")

  text(subset(labs, labs$islon)[c(2, 4), ],
       lab = parse(text = labs$lab[labs$islon][c(2, 4)]),
       pos = 1, xpd = TRUE, cex = 0.7)
  text(subset(labs, !labs$islon),
       lab = parse(text = labs$lab[!labs$islon]),
       pos = 2, xpd = TRUE, cex = 0.7)

  lines(SpPolygon)
  
  dev.off()
}

