##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Preliminary spatiotemporal plots for NBS otter taxa
##   Arctic IERP April 2023 Workshop
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(googledrive)
library(FishStatsUtils)
library(terra)
library(plotrix)
library(akgfmaps)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Spatial Objects
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbs_grid <- as.data.frame(FishStatsUtils::northern_bering_sea_grid)
nbs_pts <- terra::vect(x = nbs_grid, 
                       geom = c("Lon", "Lat"),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
nbs_pts <- terra::project(x = nbs_pts, '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

## Import land shapefile, crop to only nbs area
AK_land <- terra::vect(x = "data/spatial_data/land_shapefiles/AKland.shp")
cropped_extent <- terra::ext(nbs_pts)
ak_land_cropped <- terra::crop(x = AK_land, y = cropped_extent)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Authorize Google Drive
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
googledrive::drive_deauth()
googledrive::drive_auth()
1

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants: googledrive ids point to folder in the AIERP Synthesis 
##   Living Documents / Preliminary Analyses
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbs_otter_dir <- googledrive::as_id("1eERDrW32SGU9hABm7bRC7y5gnUAUafeg")
nbs_otter_ls <- googledrive::drive_ls(path = nbs_otter_dir)
nbs_spp <- sort(nbs_otter_ls$name)
# nbs_otter_years <- c(1985, 1988, 1991, 2010, 2017:2019, 2021:2022)
nbs_otter_years <- c(1985, 1991, 2010, 2017, 2022)
# nbs_spp <- c("Pacific cod", "walleye pollock", "yellowfin sole",
#              "bivalves", "Pacific herring",  "snow crab",
#              "Arctic cod", "saffron cod")
## Create temporary folder to put downloaded VAST files. The temp/ folder
## is in the gitignore file. 
if (!dir.exists(paths = "temp")) dir.create(path = "temp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ispp in nbs_spp) {
  
  ## Locate google drive id for ispp
  temp_spp <- nbs_otter_ls$id[nbs_otter_ls$name == ispp]
  temp_ls <- googledrive::drive_ls(path = temp_spp)
  
  ## Pull VAST fit from drive, save locally, and then load into R 
  googledrive::drive_download(file = temp_ls$id[temp_ls$name == "fit.RData"],
                              path = "temp/fit.RData",
                              overwrite = TRUE)
  
  load(file = "temp/fit.RData")
  
  ## Extract predicted densities
  D_gct <- fit$Report$D_gct[, 1, paste(nbs_otter_years)]
  
  ## Scale D_gct 
  scaled_D_gct <- D_gct / max(D_gct)
  
  ## Setup temporary shapefile with densities
  temp_nbs_pts <- nbs_pts
  temp_nbs_pts[, dimnames(D_gct)[["Time"]] ] <- scaled_D_gct
  
  jpeg(filename = paste0("presentations/2023/NPRB Arctic IERP/",
                         "April_2023_Workshop/nbs_dists/", ispp, "_nbs.jpg"), 
       width = 10, height = 1.6, units = "in", res = 500)
  par(mfrow = c(1, 6), mar = c(0, 0, 2, 0))
  for (iyear in paste0(nbs_otter_years)) {
    
    ## Rasterize the shapefile points for a given year's distribution
    nbs_ras <- terra::rast(x = temp_nbs_pts, res = 5000)
    nbs_ras <- terra::rasterize(x = temp_nbs_pts, y = nbs_ras, 
                                field = iyear)
    
    ## Plot
    image(nbs_ras, axes = F, 
          col = RColorBrewer::brewer.pal(n = 9, name = "Blues"))
    
    mtext(side = 3, text = iyear)
    plot(ak_land_cropped, add = TRUE, col = "tan", border = "tan")
  }
  
  plot(1, type = "n", axes = F, xlim = c(0, 100), ylim = c(0, 100))
  plotrix::color.legend(
    xl = 5, xr = 95,
    yb = 10, yt = 20,
    legend = round(quantile(D_gct), 0),
    rect.col = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
    gradient = "x",
    align = "rb",
    cex = 0.5)
  text(x = 50, y = 30, paste0(ispp, "\nPredicted Density (kg/km2)"), font = 2)
  dev.off()
}
