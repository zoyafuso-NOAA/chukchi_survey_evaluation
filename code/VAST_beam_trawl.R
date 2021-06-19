###############################################################################
## Project:       VAST Modelling, Chukchi Sea
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Beam Trawl Data
###############################################################################
rm(list = ls())

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
rb_data <- read.csv("data/fish_data/otter_trawl/AK_BTS_Arctic_processed.csv")
ierl_data <- read.csv("data/fish_data/2017_2019_Beam/ierl_data_processed.csv")

extrapolation_grid <- read.csv(paste0("data/spatial_data/", 
                                      "BS_Chukchi_extrapolation_grids/",
                                      "ChukchiThorsonGrid.csv"))
extrapolation_grid$Area_km2 <- extrapolation_grid$Shape_Area / 1000 / 1000

data_long <- rbind(rb_data, ierl_data[, names(rb_data)])
data_long$cpue <- data_long$catch_kg / data_long$area_swept_km2

spp_list <- c("Alaska plaice", "Arctic cod", "Bering flounder", "saffron cod",
              "snow crab", "walleye pollock" , "yellowfin sole")

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
VAST_Version <- "VAST_v12_0_0"

settings <- FishStatsUtils::make_settings( 
  Version = VAST_Version,
  n_x = 200,   # Number of knots
  Region = "User", #User inputted extrapolation grid
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits =  data.frame("STRATA" = c("All_areas"),
                              "west_border" = -Inf,
                              "east_border" = Inf), 
  bias.correct = FALSE,
  FieldConfig = c(
    "Omega1" = 1,   #Spatial random effect on occurence 
    "Epsilon1" = 1, #Spatiotemporal random effect on occurence
    "Omega2" = 1,   #Spatial random effect on positive response 
    "Epsilon2" = 1  #Spatiotemporal random effect on positive response
  ), 
  ObsModel = c(2, 1),
  max_cells = Inf,
  use_anisotropy = F,
  Options = c('SD_site_logdensity' = FALSE, 'Calculate_Range' = FALSE,
              'Calculate_effective_area' = FALSE, 'Calculate_Cov_SE' = FALSE,
              'Calculate_Synchrony' = FALSE, 'Calculate_proportion' = FALSE))


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
        "working_dir" = paste0(getwd(), "/results/beam_trawl"),
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
            file = "results/beam_trawl/model_settings.csv",
            row.names = F)
}

## Unlink Dynamic library
dyn.unload(paste0("results/beam_trawl/", VAST_Version, ".dll"))

##################################################
####   Loop over species and fit the model with "best" field configurations
##################################################

for (ispp in spp_list[-1]) { ## Loop over species -- start
  
  ## Subset data input for species ispp
  data_geostat_subset <- subset(x = data_geostat, 
                                subset = spp == ispp)
  
  ## Subset the model runs that converged with a low enough maximum gradient 
  sub_df <- subset(x = model_settings, 
                   subset = common_name == ispp & max_grad < 1e-4)
  
  ## If there is at least one converged model, run the model with the lowest
  ## RRMSE in the density predictions
  if (nrow(sub_df) > 0) {
    field_config <- unlist(sub_df[which.min(sub_df$rrmse), 
                                               c("Omega1", "Epsilon1", 
                                                 "Omega2", "Epsilon2")])
    settings$FieldConfig <- field_config
    
    ## Fit model
    result_dir <- paste0(getwd(), "/results/beam_trawl/", ispp, "/")
    fit <- FishStatsUtils::fit_model( 
      "settings" = settings,
      "working_dir" = result_dir,
      "Lat_i" = data_geostat_subset[, "Lat"],
      "Lon_i" = data_geostat_subset[, "Lon"],
      "t_i" = data_geostat_subset[, "Year"],
      "b_i" = data_geostat_subset[, "Catch_KG"],
      "a_i" = data_geostat_subset[, "AreaSwept_km2"],
      "getJointPrecision" = TRUE,
      "newtonsteps" = 1,
      "test_fit" = F,
      "input_grid" = extrapolation_grid)
    
    ## Diagnostics
    diagnostics <- plot(fit, working_dir = result_dir)
    
    ## Save fit and diagnostics
    fit_sub <- fit$Report
    save(list = "fit_sub", file = paste0(result_dir, "/fit.RData"))
    save(list = "diagnostics", file = paste0(result_dir, "/diagnostics.RData"))
    
    ## Unlink Dynamic library
    dyn.unload(paste0(result_dir, "/", VAST_Version, ".dll"))
  }
  
} ## Loop over species -- end
