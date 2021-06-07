###############################################################################
## Project:       VAST Modelling, Chukchi Sea
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   2012 Chukchi data
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
wd <- "C:/Users/zack.oyafuso/Desktop/Arctic/"
output_wd <- paste0(wd, "VAST_model_output_2012_2019_BT/")

##################################################
####  Import Libraries
##################################################
library(tidyr)
library(reshape)
library(VAST)
library(raster)
library(sp)
library(readxl)
library(rgdal)

##################################################
####  Import Data
##################################################
rb_data <- read.csv(paste0(wd, "data/AK_BTS_ARCTIC.csv"))
ierl_data <- read.csv(paste0(wd, "data/ierl_data.csv"))

extrapolation_grid <- read.csv(paste0(wd, "spatial_data/", 
                                      "BS_Chukchi_extrapolation_grids/",
                                      "ChukchiThorsonGrid.csv"))
extrapolation_grid$Area_km2 <- extrapolation_grid$Shape_Area / 1000 / 1000

data_long <- rbind(rb_data, ierl_data[, names(rb_data)])
data_long$cpue <- data_long$catch_kg / data_long$area_swept_km2

# no_recs <- spread(data = aggregate(cpue ~  year + common_name,
#                                    data = data_long,
#                                    subset = gear == "beam",
#                                    FUN = function(x) sum(x > 0, na.rm = T)),
#                   value = "cpue",
#                   key = "year"
# )
# 
# spp_no_recs <- no_recs$common_name[apply(X = no_recs[, c("2012",
#                                                          "2017",
#                                                          "2019")],
#                                          MARGIN = 1,
#                                          FUN = function(x) !all(x > 10,
#                                                                 na.rm = TRUE))]
# 
# spp_list <- no_recs$common_name[!no_recs$common_name %in% spp_no_recs]

spp_list <- c(#"Alaska plaice",
              "Arctic cod", "Bering flounder", "saffron cod", "snow crab", 
              "walleye pollock" , "yellowfin sole",
              "Arctic staghorn sculpin", "circumboreal toad crab",
              "fuzzy hermit crab", "hairy hermit crab", "northern nutclam",
              "notched brittlestar", "shorthorn (=warty) sculpin",
              "slender eelblenny"
              )

data_long <- subset(x = data_long,
                    subset = gear == "beam" & 
                      common_name %in% spp_list &
                      year %in% c(2012, 2017, 2019))

data_geostat <- data.frame( 
  spp = as.factor(data_long$common_name),
  Year = data_long$year,
  Catch_KG = data_long$cpue,
  AreaSwept_km2 = 1,
  Lat = data_long$lat,
  Lon = data_long$lon, 
  stringsAsFactors = T)

##################################################
####  Model settings
##################################################
settings <- FishStatsUtils::make_settings( 
  n_x = 100,   # Number of knots
  Region = "User", #User inputted extrapolation grid
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits =  data.frame("STRATA" = c("All_areas"),
                              "west_border" = -Inf,
                              "east_border" = Inf), 
  bias.correct = FALSE,
  FieldConfig = c(
    "Omega1" = 3,   #Spatial random effect on occurence 
    "Epsilon1" = 3, #Spatiotemporal random effect on occurence
    "Omega2" = 3,   #Spatial random effect on positive response 
    "Epsilon2" = 3  #Spatiotemporal random effect on positive response
  ), 
  RhoConfig = c("Beta1" = 0, 
                "Beta2" = 0, 
                "Epsilon1" = 0, 
                "Epsilon2" = 0), #Each year is a fixed effect
  OverdispersionConfig = c("Eta1" = 0, 
                           "Eta2" = 0), #Turn off overdispersion 
  "Options" = c("Calculate_Range" = F, 
                "Calculate_effective_area" = F),
  ObsModel = c(2, 1),
  max_cells = Inf,
  use_anisotropy = F)


##################################################
####   Model result objects
##################################################
model_settings <- expand.grid(Omega1 = 1,
                              Omega2 = 1,
                              Epsilon1 = 0:1,
                              Epsilon2 = 0:1,
                              common_name = spp_list,
                              stringsAsFactors = FALSE)

model_settings[, c("status", "max_grad", "rrmse")] <- NA

##################################################
####   Model fit
##################################################
for (irow in 1:nrow(model_settings)) {
  settings$FieldConfig <- unlist(model_settings[irow, c("Omega1", "Epsilon1",
                                                        "Omega2", "Epsilon2")])
  data_geostat_subset <- subset(x = data_geostat, 
                                subset = spp == model_settings$common_name[irow])
  
  fit <- tryCatch(
    {
      FishStatsUtils::fit_model( 
        "settings" = settings,
        "working_dir" = output_wd,
        "Lat_i" = data_geostat_subset[, "Lat"],
        "Lon_i" = data_geostat_subset[, "Lon"],
        "t_i" = data_geostat_subset[, "Year"],
        "b_i" = data_geostat_subset[, "Catch_KG"],
        "a_i" = data_geostat_subset[, "AreaSwept_km2"],
        "getJointPrecision" = TRUE,
        "newtonsteps" = 1,
        "test_fit" = F,
        "input_grid" = extrapolation_grid)
    },
    error = function(cond) {
      message("Did not converge. Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message( paste("Processed", model_settings$common_name[irow],
                     "with settings",
                     "Omega1 =", model_settings$Omega1[irow],
                     "Omega2 =", model_settings$Omega2[irow],
                     "Epsilon1 =", model_settings$Epsilon1[irow],
                     "Epsilon2 =", model_settings$Epsilon2[irow]) )
    }
  )    
  
  ##################################################
  ####   Did the model fit fun produce an error?
  ##################################################
  model_settings$status[irow] <- ifelse(test = is.na(fit), 
                                        yes = "no_convergence",
                                        no = "confirm_gradient")
  
  if (model_settings$status[irow]  == "confirm_gradient") {
    
    ## Maximum Gradient
    model_settings$max_grad[irow] <- fit$parameter_estimates$max_gradient
    
    ## RRMSE of density predictions
    obs_cpue <- with(data_geostat_subset,  Catch_KG / AreaSwept_km2 )
    pred_cpue <- fit$Report$D_i
    rrmse <- sqrt(mean((obs_cpue - pred_cpue)^2)) / mean(obs_cpue)
    
    model_settings$rrmse[irow] <- round(rrmse, 3)
  } else (model_settings[irow, c("status", "rrmse", "max_grad")] <- NA)
  
  ##################################################
  ####   Save
  ##################################################
  write.csv(model_settings, 
            file = paste0(output_wd, "model_settings.csv"),
            row.names = F)
}


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

cov2cor(fit$Report$L_epsilon2_cf %*% t(fit$Report$L_epsilon2_cf))
cov2cor(fit$Report$L_epsilon1_cf %*% t(fit$Report$L_epsilon1_cf))

cov2cor(fit$Report$L_omega1_cf %*% t(fit$Report$L_omega1_cf))
cov2cor(fit$Report$L_omega2_cf %*% t(fit$Report$L_omega2_cf))



