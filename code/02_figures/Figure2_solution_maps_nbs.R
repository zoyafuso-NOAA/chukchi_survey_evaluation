##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Figure X: Plot STRS optimized solutions
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot STRS optimization solutions with different allocations
##                     based on species tradeoffs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import relevant libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)
library(graticule)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##   Crop the Alaska shapefile to only include NW AK
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
load(here::here("data/survey_opt_data/optimization_data.RData"))
ak_land <- rgdal::readOGR(here::here("data/spatial_data/land_shapefiles/AKland.shp"))
cropped_extent <- extent(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- raster::crop(x = ak_land, y = cropped_extent)

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
##   Plot ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
{
  ## Open png device
  png(filename = here::here(paste0("figures/Figure2_ms_strs_sols_nbs.png")),
      units = "mm", width = 90, height = 125, res = 500)
  
  ## Panel Layout
  par(mfrow = c(2, 2),
      mar = c(0, 0, 0, 0), 
      oma = c(0.5, 0.5, 0.5, 0.5), 
      family = "serif")
  
  for (igear in c("otter", "beam")) { ## Loop over gears -- start
    for (istrata in 3:4) { ## Loop over strata -- start
      
      load(here::here(paste0("results/nbs_", igear, "/survey_opt/Str_", 
                  istrata, "/result_list_nbs.RData")))
      load(here::here(paste0("results/nbs_", igear, "/survey_opt/Str_",
                  istrata, "/allocations_nbs.RData")))
      
      ## solution shape object
      solution <- sp::SpatialPointsDataFrame(
        coords = grid_pts@coords,
        data = data.frame(Str_no = result_list$sol_by_cell), 
        proj4string = aea_crs )
      solution_ras <- raster::raster(x = solution,
                                     resolution = 4000)
      solution_ras <- raster::rasterize(x = solution,
                                        y = solution_ras,
                                        field = "Str_no")
      strata_cols <- RColorBrewer::brewer.pal(n = istrata, name = "Set2") 
      
      ## 100-sample allocations
      ms_solutions <- subset(x = ms_sample_allocations,
                             subset = n == 100,
                             select = paste0("Str_", 1:istrata))
      
      ## Random draw from the multispecies allocation
      sample_vec <- c()
      for(i in 1:istrata) {
        sample_vec <- c(sample_vec,
                        sample(x = which(result_list$sol_by_cell == i),
                               size = ms_solutions[, i]) )
      }
      
      ## Turn randomly drawn stations into a spatial object
      samp_pts <- sp::SpatialPoints(
        coords = grid_pts@coords[sample_vec, ], 
        proj4string = aea_crs)
      
      ## Plot solution
      image(x = solution_ras, 
            col = strata_cols, asp = 1, axes = F, ann = F)
      
      ## Plot sampled stations
      points(samp_pts, pch = 16, cex = 0.75)
      
      ## Plot graticules
      plot(aea_grat, add = T, lwd = 0.5)
      
      ## Plot land
      plot(ak_land_cropped, col = "tan", border = F, add = TRUE)
      
      ## Plot graticule labels
      text(subset(aea_grat_labs, aea_grat_labs$islon), 
           lab = parse(text = aea_grat_labs$lab[aea_grat_labs$islon]), 
           pos = 3, cex = 0.7)
      text(subset(aea_grat_labs, !aea_grat_labs$islon), 
           lab = parse(text = aea_grat_labs$lab[!aea_grat_labs$islon]),
           pos = 1, cex = 0.7)
      
      ## Plot gear/stratum label, sampling proportions across strata
      legend(x = extent(x = solution_ras)[1] + 
               diff(x = extent(x = solution_ras)[1:2]) * 0.5,
             y = extent(x = solution_ras)[3] + 
               diff(x = extent(x = solution_ras)[3:4]) * 0.5,
             legend = (x = paste0(
               "Str ", 1:istrata, ": ", 
               ms_solutions[, paste0("Str_", istrata:1)]/100)),
             title = paste0(igear, " trawl"),
             fill = rev(strata_cols), cex = 0.85, bty = "n")
      
      ## Scale bar
      segments(x0 = extent(solution_ras)[1] + 
                 diff(extent(solution_ras)[1:2])*0.75, 
               x1 = extent(solution_ras)[1] + 
                 diff(extent(solution_ras)[1:2])*0.75 + 100000,
               y0 = extent(solution_ras)[3] + 
                 diff(extent(solution_ras)[3:4])*0.1,
               lwd = 2)
      text(x = extent(solution_ras)[1] + 
             diff(extent(solution_ras)[1:2]) * 0.83,
           y = extent(solution_ras)[3] + 
             diff(extent(solution_ras)[3:4]) * 0.1,
           labels = "100 km", pos = 1,
           cex = 0.9)
      
      box()
      
    } ## Loop over strata -- start
  } ## Loop over gears -- end
  
  ## Close device
  dev.off()
}
