##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Figure 1: Plot STRS optimized solutions
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot STRS optimization solutions with different allocations
##                     based on species tradeoffs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import relevant libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(terra)
library(plotrix)
library(RColorBrewer)
library(graticule)
library(FishStatsUtils)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##   Crop the Alaska shapefile to only include NW AK
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spp_settings <- read.csv(file = here::here("results/good_species.csv"), 
                         stringsAsFactors = F)
load(file = here::here("data/survey_opt_data/optimization_data.RData"))
ak_land <- terra::vect(x = here::here("data/spatial_data/land_shapefiles/AKland.shp"))
cropped_extent <- terra::ext(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- terra::crop(x = ak_land, y = cropped_extent)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set up graticules ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
aea_grat <- graticule::graticule(
  lons = seq(from = -180, to = -150, by = 5),
  lats = seq(from = 65, to = 75, by = 2.5), 
  proj = aea_crs)
aea_grat_labs <- graticule::graticule_labels(
  lons = seq(from = -180, to = -150, by = 10),
  lats = seq(from = 70, to = 72.5, by = 2.5),
  xline = -171, yline = 72.5,
  proj = aea_crs)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Which species were represented by both gears or just one gear?
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
col_otter <- c("white", brewer.pal(n = 9, name = "Blues"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import fitted densities for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D_gct <- array(dim = c(n_cells, n_spp_beam + n_spp_otter),
               dimnames = list(NULL, 
                               c(paste0("beam_", spp_list_beam), 
                                 paste0("otter_", spp_list_otter))))

for (igear in c("beam", "otter")) { ## Loop over gears -- start
  temp_spp_list <- spp_list$taxon
  
  ## index for the most recent year of data
  iyear <- switch(igear, "beam" = 8, "otter" = length(1985:2022))

  for (ispp in temp_spp_list) { ## Loop over species -- start
    
    ## Load VAST fit for the gear/taxon combo
    load(here::here(paste0("results/nbs_", igear, "/vast_fits/", ispp, "/fit.RData")))
    
    ## Record the predicted density for the most recent survey year
    D_gct[, paste0(igear, "_", ispp)] <- fit$Report$D_gct[, 1, iyear]
    
  } ## Loop over species -- end
} ## Loop over gears -- end

dimnames(D_gct)

## append locations to D_gct
D_gct_w_pts <- data.frame(northern_bering_sea_grid,
                          data.frame(D_gct, check.names = F),
                          check.names = F)

## Project to aea by projecting from lon/lat
D_gct_pts <- terra::vect(x = D_gct_w_pts,
                         geom = c("Lon", "Lat"),
                         crs = latlon_crs)
D_gct_pts <- terra::project(x = D_gct_pts, aea_crs)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot ----
##   There are 18 species just represented by the otter trawl 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
{
  # Open png device
  png(filename = paste0("figures/Figure1_pred_density.png"),
      units = "mm", width = 190, height = 200, res = 500)
  
  ## Set the layout of the plot
  layout(mat = matrix(data = c(1, 2, 30, 3, 4, 30, 5, 6, 30,
                               7, 8, 30, 9,10, 30, 11,12,30, 
                               13,14,30, 15,16,30, 17,18,30, 
                               19,20,30, 21,22,30, 23,24,30,
                               26,27,30, 28,29,30, 30,25,30),
                      nrow = 5, byrow = TRUE),
         widths = c(1,1,0.2, 1,1,0.2, 1,1,0.2))
  
  ## Set plot width margins
  par(mar = c(0, 0, 2.5, 0),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = "serif")
  
  for (ispp in c(otter_beam_idx, 
                 otter_spp_idx, 
                 beam_spp_idx)) { ## Loop over species -- start
    for (igear in c("otter", "beam")) { ## Loop over gears -- start
      
      if (igear == "otter" & ispp %in% beam_spp_idx) next
      if (igear == "beam" & ispp %in% otter_spp_idx) next
      
      ## Create solution shape object
      temp_raster <- terra::rast(terra::ext(x = D_gct_pts), 
                                 resolution = c(5000, 5000) )
      
      temp_raster <- terra::rasterize(x = D_gct_pts, 
                                      y = temp_raster,
                                      field = paste0(igear, "_", 
                                                     spp_settings$taxon[ispp]))
      
      ## Plot solution
      image(x = temp_raster, asp = 1, axes = F, ann = F,
            col = get(paste0("col_", igear)))
      
      ## Plot graticules
      plot(aea_grat, add = T, lwd = 0.5)
      
      ## Plot land
      plot(ak_land_cropped, col = "tan", border = F, add = TRUE)
      
      ## Plot species title
      if (ispp %in% otter_beam_idx & igear == "otter")
        text(x = ext(temp_raster)[1] +
               diff(ext(temp_raster)[1:2]) * 1,
             y = ext(temp_raster)[3] +
               diff(ext(temp_raster)[3:4]) * 1.225,
             labels = gsub(pattern = '\\n',
                           x = spp_settings$taxon_plot[ispp],
                           replacement = "\n",
                           fixed = TRUE),
             pos = 1, xpd = NA, cex = 1.5, font = 2)
      
      if (ispp %in% c(otter_spp_idx, beam_spp_idx))
        mtext(side = 3, line = 0.75,
              text = spp_settings$taxon_plot[ispp],
              cex = 1, font = 2)
      
      ## Plot graticule labels
      text(subset(x = aea_grat_labs, 
                  subset = lab %in% paste0(c(170,160), "*degree*W")),
           lab = parse(text = aea_grat_labs$lab[aea_grat_labs$lab %in% 
                                                  paste0(c(170, 160), 
                                                         "*degree*W")]),
           pos = 3, cex = 0.65, xpd = NA)
      text(subset(aea_grat_labs, !aea_grat_labs$islon), 
           lab = parse(text = aea_grat_labs$lab[!aea_grat_labs$islon]),
           pos = 1, cex = 0.65, xpd = NA)
      
      ## Scale bar
      segments(x0 = ext(x = temp_raster)[1] + 
                 diff(x = ext(x = temp_raster)[1:2]) * 0.75, 
               x1 = ext(x = temp_raster)[1] + 
                 diff(x = ext(x = temp_raster)[1:2]) * 0.75 + 100000,
               y0 = ext(x = temp_raster)[3] + 
                 diff(x = ext(x = temp_raster)[3:4]) * 0.1,
               lwd = 2)
      text(x = ext(x = temp_raster)[1] + 
             diff(x = ext(x = temp_raster)[1:2]) * 0.83,
           y = ext(x = temp_raster)[3] + 
             diff(x = ext(x = temp_raster)[3:4]) * 0.11,
           labels = "100 km", 
           pos = 1, cex = 0.75)
      
      ## Density Legend
      zlim_ <- as.numeric(x = terra::minmax(x = temp_raster))
      
      plotrix::color.legend(
        yb = ext(x = temp_raster)[3] + 
          diff(x = ext(x = temp_raster)[3:4]) * 0.15,
        yt = ext(x = temp_raster)[3] + 
          diff(x = ext(x = temp_raster)[3:4]) * 0.55,
        xl = ext(x = temp_raster)[1] + 
          diff(x = ext(x = temp_raster)[1:2]) * 0.75,
        xr = ext(x = temp_raster)[1] + 
          diff(x = ext(x = temp_raster)[1:2]) * 0.8,
        legend = pretty(zlim_, n = 3),
        rect.col = get(paste0("col_", igear)),
        gradient = "y", 
        align = "rb",
        cex = 0.5)
      
      ## Gear Label
      text(x = ext(x = temp_raster)[1] + 
             diff(x = ext(x = temp_raster)[1:2]) * 0.56,
           y = ext(x = temp_raster)[3] + 
             diff(x = ext(x = temp_raster)[3:4]) * 0.4,
           labels = paste0(igear, " trawl\n(kg/km2)"), 
           pos = 1, cex = 0.75)
      
    } ## Loop over gears -- end
  } ## Loop over species -- end
  dev.off()
}
