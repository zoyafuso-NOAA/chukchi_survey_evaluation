
##################################################
####   Mapping function
##################################################
make_a_raster <- function(extrap_grid = extrapolation_grid,
                          lat_lon_names = c("Lon", "Lat"),
                          plot_what = fit$Report$D_gct[, 1, 1],
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

###############################################################
par(mfcol = c(2, length(spp_list)), 
    mar = c(2.5, 3.5, 2.5, 3.5), 
    oma = c(0, 0, 1, 0))
for(iyear in c(1, 3)) {
  for (which_spp_code in 1:length(spp_list)) {
    
    which_spp <- sort(spp_list)[which_spp_code]
    plot_this <- fit$Report$D_gct[, which_spp_code, iyear]
    
    density_ras <- make_a_raster(plot_what = plot_this)
    plot(density_ras, 
         las = 1, 
         main = ifelse(iyear == 1, which_spp, NA),
         cex.main = 2)
    
    temp_df <- subset(x = data_wide,
                      subset = year == c(2017:2019)[iyear],
                      select = c("lat", "lon", "year", 
                                 which_spp))
    points(lat ~ lon,
           data = temp_df[temp_df[, which_spp] == 0, ],
           cex = 1,
           pch = 16,
           col = 'red')
    points(lat ~ lon,
           data = temp_df,
           pch = 1,
           cex = temp_df[, which_spp] / max(temp_df[, which_spp]) * 5 )
    # plot(AK,
    #      add = TRUE,
    #      col = "tan",
    #      border = F)
  }
}

par(mfcol = c(2, length(spp_list)), mar = c(3, 3.5, 2, 3.5))
for(ipred in 1:2) {
  for (which_spp_code in 1:length(spp_list)  ) {
    which_spp <- sort(spp_list)[which_spp_code]
    plot_this <- fit$Report[[paste0("Omega", 
                                    ipred, 
                                    "_gc")]][, which_spp_code]
    
    density_ras <- make_a_raster(plot_what = plot_this)
    
    plot(density_ras, 
         las = 1, 
         main = ifelse(ipred == 1, which_spp, NA),
         cex.main = 2)
    plot(AK,
         add = TRUE,
         col = "tan",
         border = F)
  }
}

par(mfcol = c(4, length(spp_list)), mar = c(3, 3.5, 2, 3.5))
for(iyear in c(1, 3)) {
  for(ipred in 1:2) {
    for (which_spp_code in 1:length(spp_list)  ) {
      which_spp <- sort(spp_list)[which_spp_code]
      plot_this <- fit$Report[[paste0("Epsilon", 
                                      ipred, 
                                      "_gct")]][, which_spp_code, iyear]
      
      density_ras <- make_a_raster(plot_what = plot_this)
      
      plot(density_ras, 
           las = 1, 
           main = ifelse(ipred == 1, which_spp, NA),
           cex.main = 2)
      # plot(AK,
      #      add = TRUE,
      #      col = "tan",
      #      border = F)
    }
  }
}
