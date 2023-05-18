##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Figure X: Plot STRS optimized solutions
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot STRS optimization solutions with different allocations
##                     based on species tradeoffs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import relevant libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(terra)
library(RColorBrewer)
library(graticule)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##   Crop the Alaska shapefile to only include NW AK
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
load(file = "data/survey_opt_data/optimization_data.RData")
ak_land <- terra::vect(x = "data/spatial_data/land_shapefiles/AKland.shp")
grid_pts <- terra::vect(x = "data/survey_opt_data/grid_pts.shp")
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
##   Plot ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
{
  ## Open jpeg device
  jpeg(filename = paste0("figures/Figure2_ms_strs_sols.jpeg"),
      units = "mm", width = 90, height = 125, res = 500)
  
  ## Panel Layout
  par(mfrow = c(2, 2),
      mar = c(0, 0, 0, 0), 
      oma = c(0.5, 0.5, 0.5, 0.5), 
      family = "serif")
  
  ipanel <- 1
  for (igear in c("otter", "beam")) { ## Loop over gears -- start
    for (istrata in 3:4) { ## Loop over strata -- start
      
      load(file = paste0("results/chukchi_", igear, "/survey_opt/Str_", 
                         istrata, "/result_list.RData"))
      load(file = paste0("results/chukchi_", igear, "/survey_opt/Str_",
                         istrata, "/allocations.RData"))
      
      ## solution shape object
      solution <- 
        terra::vect(x = data.frame(cbind(terra::crds(grid_pts),
                                         Str_no = result_list$sol_by_cell)),
                    geom = c("x", "y"),
                    crs = aea_crs)
      solution_ras <- terra::rast(x = solution, resolution = 4000)
      solution_ras <- terra::rasterize(x = solution,
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
      samp_pts <- terra::vect(x = terra::crds(grid_pts)[sample_vec, ],
                              crs = aea_crs)

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
      legend(x = terra::ext(x = solution_ras)[1] + 
               diff(x = terra::ext(x = solution_ras)[1:2]) * 0.5,
             y = terra::ext(x = solution_ras)[3] + 
               diff(x = terra::ext(x = solution_ras)[3:4]) * 0.5,
             legend = (x = paste0(
               "Str ", 1:istrata, ": ", 
               ms_solutions[, paste0("Str_", istrata:1)]/100)),
             title = paste0(LETTERS[ipanel], ") ", igear, " trawl"),
             fill = rev(strata_cols), cex = 0.85, bty = "n")
      
      ## Scale bar
      segments(x0 = terra::ext(solution_ras)[1] + 
                 diff(terra::ext(solution_ras)[1:2])*0.75, 
               x1 = terra::ext(solution_ras)[1] + 
                 diff(terra::ext(solution_ras)[1:2])*0.75 + 100000,
               y0 = terra::ext(solution_ras)[3] + 
                 diff(terra::ext(solution_ras)[3:4])*0.1,
               lwd = 2)
      text(x = terra::ext(solution_ras)[1] + 
             diff(terra::ext(solution_ras)[1:2]) * 0.83,
           y = terra::ext(solution_ras)[3] + 
             diff(terra::ext(solution_ras)[3:4]) * 0.1,
           labels = "100 km", pos = 1,
           cex = 0.9)
      
      box()
      
      ipanel <- ipanel + 1
    } ## Loop over strata -- start
  } ## Loop over gears -- end
  
  ## Close device
  dev.off()
}
