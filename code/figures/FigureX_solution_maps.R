###############################################################################
## Project:       Figure X: Plot STRS optimized solutions
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot STRS optimization solutions with different allocations
##                     based on species tradeoffs
###############################################################################
rm(list = ls())

##################################################
#### Import relevant libraries   
##################################################
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)
library(graticule)

##################################################
#### Import Data and Constants
#### Crop the Alaska shapefile to only include NW AK
##################################################
load("data/survey_opt_data/optimization_data.RData")
ak_land <- rgdal::readOGR("data/spatial_data/land_shapefiles/AKland.shp")
cropped_extent <- extent(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- raster::crop(x = ak_land, y = cropped_extent)

##################################################
#### Graticules
##################################################
aea_grat <- graticule::graticule(lons = seq(from = -175, to = -155, by = 5),
                                 lats = seq(from = 65, to = 75, by = 2.5), 
                                 proj = aea_crs)
aea_grat_labs <- graticule::graticule_labels(lons = seq(from = -175, 
                                                        to = -160, by = 5),
                                             lats = seq(from = 70, 
                                                        to = 72.5, by = 2.5),
                                             xline = -171, yline = 72.5,
                                             proj = aea_crs)

##################################################
#### 
##################################################

for (igear in c("otter", "beam")) { ## Loop over gears -- start
  
  ## Open png device
  png(filename = paste0("figures/FigureX_", igear, "_sols.png"),
      units = "mm", width = 190, height = 100, res = 100)
  
  ## Set species and order depending on the gear
  spp_list <- get(paste0("spp_list_", igear))
  ns <- length(spp_list)
  spp_order <- list("otter" = c(2, 5, 3, 4, 6, 1),
                    "beam" = c(1, 4, 2, 3, 5))[[igear]]
  panel_order <- list("otter" = c(1, 2, 3, 7, 7, 
                                  4, 5, 6, 7, 7),
                      "beam" = c(1, 2, 3, 6, 6, 
                                 4, 5, 7, 6, 6))[[igear]]
  
  ## Load 4-strata data
  istrata = 4
  load(paste0("results/chukchi_", igear, "/survey_opt/MS/Str_", 
              istrata, "/result_list.RData"))
  load(paste0("results/chukchi_", igear, "/survey_opt/MS/Str_", 
              istrata, "/allocations.RData"))
  
  ## solution shape object
  solution <- sp::SpatialPointsDataFrame(
    coords = grid_pts@coords,
    data = data.frame(Str_no = result_list$sol_by_cell), 
    proj4string = aea_crs )
  solution_ras <- raster::raster(x = solution,
                                 resolution = 3700)
  solution_ras <- raster::rasterize(x = solution,
                                    y = solution_ras,
                                    field = "Str_no")
  strata_cols <- RColorBrewer::brewer.pal(n = istrata, name = "Set2") 
  
  ## 100-sample allocations
  ss_solutions <- subset(x = ss_sample_allocations, 
                         subset = n == 100, 
                         select = paste("Str_", 1:istrata))
  ms_solutions <- subset(x = ms_sample_allocations,
                         subset = n == 100,
                         select = paste("Str_", 1:istrata))
  
  ## Panel Layout
  layout(mat = matrix(data = panel_order,
                      byrow = TRUE, nrow = 2))
  par(mar = c(0, 0, 0, 0), oma = c(0.5, 0.5, 0.5, 0.5), family = "serif")
  
  ## Loop over species and plot solution and sample allocation
  for (ispp in spp_order) { ## Loop over species -- start
    
    ## Sample from the allocation
    sample_vec <- c()
    for(i in 1:istrata) {
      sample_vec <- c(sample_vec,
                      sample(x = which(result_list$sol_by_cell == i),
                             size = ss_solutions[ispp, i]) )
    }
    
    ## Turn the sampled stations into a spatial object
    samp_pts <- sp::SpatialPoints(
      coords = grid_pts@coords[sample_vec, ], 
      proj4string = aea_crs)
    
    ## Plot solution
    image(x = solution_ras, 
          col = strata_cols,
          asp = 1,
          axes = F, ann = F)
    
    ## Plot sampled stations
    points(samp_pts,
           pch = 16, cex = 0.5)
    
    ## Plot graticules
    plot(aea_grat, add = T)
    
    ## Plot land
    plot(ak_land_cropped, col = "tan", border = F, add = TRUE)
    
    ## Plot species label
    text(x = extent(solution_ras)[1] + diff(extent(solution_ras)[1:2])*0.75,
         y = extent(solution_ras)[3] + diff(extent(solution_ras)[3:4])*0.3,
         labels = gsub(x = spp_list[ispp], pattern = " ", replacement = "\n"),
         cex = 1.5)
    box()
  }  ## Loop over species -- end
  
  ## Multispecies allocation
  sample_vec <- c()
  for(i in 1:istrata) {
    sample_vec <- c(sample_vec,
                    sample(x = which(result_list$sol_by_cell == i),
                           size = ms_solutions[, i]) )
  }
  
  samp_pts <- sp::SpatialPoints(
    coords = grid_pts@coords[sample_vec, ], 
    proj4string = aea_crs)
  
  ## Plot solution
  image(x = solution_ras, 
        col = strata_cols,
        asp = 1,
        axes = F, ann = F)
  
  ## Plot sampled stations
  points(samp_pts,
         pch = 16, cex = 0.75)
  
  ## Plot graticules
  plot(aea_grat, add = T)
  
  ## Plot land
  plot(ak_land_cropped, col = "tan", border = F, add = TRUE)
  
  ## Plot species label
  text(x = extent(solution_ras)[1] + diff(extent(solution_ras)[1:2])*0.75,
       y = extent(solution_ras)[3] + diff(extent(solution_ras)[3:4])*0.3,
       labels = "Multispecies\nSolution", cex = 2 )
  
  ## Plot graticule labels
  text(subset(aea_grat_labs, aea_grat_labs$islon), 
       lab = parse(text = aea_grat_labs$lab[aea_grat_labs$islon]), pos = 3)
  text(subset(aea_grat_labs, !aea_grat_labs$islon), 
       lab = parse(text = aea_grat_labs$lab[!aea_grat_labs$islon]),
       pos = 1)
  box()
  
  ## Close device
  dev.off()
} ## Loop over gears -- end
