###############################################################################
## Project:       Plot function 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create function to output basic solution maps
###############################################################################

plot_survey_opt_map <- function( file_name,
                                 grid_object,
                                 sol_by_cell,
                                 draw_stations,
                                 allocations){
  
  ## Import libraries
  library(sp); library(raster)
  
  ## Set up spatial object
  goa <- sp::SpatialPointsDataFrame(
    coords = grid_object@coords,
    data = data.frame(Str_no = sol_by_cell), 
    proj4string = crs(grid_object) )
  goa_ras <- raster::raster(x = goa,
                            resolution = 5000)
  goa_ras <- raster::rasterize(x = goa,
                               y = goa_ras,
                               field = "Str_no")
  
  ## Set up plot
  png(filename = file_name, width = 6, height = 3, units = "in", res = 500)
  
  ## Set up panel layout
  par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
  
  ## Plot solution
  strata_colors <- colorRampPalette(
    RColorBrewer::brewer.pal(n = 11,
                             name = "Paired"))(length(unique(plot_solution)) )
  
  if (any(sol_by_cell == 0)) strata_colors <- c("black", strata_colors)
  
  plot(goa_ras, axes = F, asp = 1, col =  strata_colors)
  
  ## Draw samples if needed
  if (draw_stations) {
    
    ## Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:length(allocations)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(sol_by_cell == istrata),
                             size = allocations[istrata]) )
    }
    
    samp_pts <- sp::SpatialPoints(
      coords = grid_object@coords[sample_vec, ], 
      proj4string = crs(grid_object))
    
    points(samp_pts,
           pch = 16, cex = 0.5)
  }
  
  ## Close device
  dev.off()
}
