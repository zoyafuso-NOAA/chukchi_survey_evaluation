###############################################################################
## Project:       VAST Modelling, Chukchi Sea
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   1990 and 2012 Chukchi data
## 
## Notes:         R Version 4.0.2
##                VAST 3.6.1; FishStatsUtils 2.8.10
###############################################################################
rm(list = ls())
r_version <- "R version 4.0.2 (2020-06-22)"
vast_version <- "3.6.1"
fishstatsutils_version <- "2.8.0"
vast_cpp_version <- "VAST_v12_0_0"

ifelse(test = sessionInfo()$R.version$version.string == r_version &
         packageVersion("VAST") == vast_version &
         packageVersion("FishStatsUtils") == fishstatsutils_version &
         FishStatsUtils::get_latest_version() == vast_cpp_version,
       yes = "versions are good to go", 
       no = "check your versions")


##################################################
####  Import Libraries
##################################################
# library(tidyr)
# library(reshape)
# 
# library(raster)
# library(sp)
# library(rgdal)

# library(devtools)
# devtools::install_local("C:/Users/zack.oyafuso/Downloads/FishStatsUtils-2.8.0")
# devtools::install_github("James-Thorson-NOAA/VAST@3.6.1") #for a specific VAST package version
library(VAST)

##################################################
####  Import CPUE Data and Grid
##################################################
spp_list <- c("Alaska plaice", "Arctic cod", "Bering flounder", "saffron cod",
              "walleye pollock", "yellowfin sole", "snow crab")

rb_data <- read.csv("data/fish_data/otter_trawl/AK_BTS_Arctic_processed.csv")

extrapolation_grid <- read.csv(paste0("data/spatial_data/",
                                      "BS_Chukchi_extrapolation_grids/",
                                      "ChukchiThorsonGrid.csv"))
extrapolation_grid$Area_km2 <- extrapolation_grid$Shape_Area / 1000 / 1000

data_long <- subset(x = rb_data, 
                    subset = gear == "otter" & common_name %in% spp_list) 
data_long$cpue <- data_long$catch_kg / data_long$area_swept_km2
data_long <- na.omit(data_long)

data_geostat <- data.frame( 
  spp = as.factor(data_long$common_name),
  Year = data_long$year,
  Catch_KG = data_long$cpue,
  AreaSwept_km2 = 1,
  Lat = data_long$lat,
  Lon = data_long$lon,
  Depth = scale(data_long$bot_depth),
  Temp = scale(data_long$bot_temp),
  stringsAsFactors = T)

##################################################
####   VAST Model Settings
##################################################
VAST_Version <- "VAST_v12_0_0"
settings <- FishStatsUtils::make_settings( 
  Version = VAST_Version,
  n_x = 200,   # Number of knots
  Region = "User", #User inputted extrapolation grid
  purpose = "index2",
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
model_settings <- expand.grid(common_name = spp_list,
                              Omega1 = 0:1,
                              Omega2 = 0:1,
                              Epsilon1 = 0:1,
                              Epsilon2 = 0:1,
                              dens_covar = c("FALSE"),
                              stringsAsFactors = FALSE)

model_settings[, c("status", "max_grad", "rrmse", "aic")] <- NA

##################################################
####   Create result directory if not created already
##################################################
if(!dir.exists("results/otter_trawl/")) dir.create("results/otter_trawl/")

##################################################
####   Model fit
##################################################
for (irow in 1:nrow(model_settings)) {
  settings$FieldConfig <- unlist(model_settings[irow, c("Omega1", "Epsilon1",
                                                        "Omega2", "Epsilon2")])
  
  data_geostat_subset <- subset(x = data_geostat, 
                                subset = spp == model_settings$common_name[irow])
  
  fit <- tryCatch( {FishStatsUtils::fit_model( 
    "settings" = settings,
    "working_dir" = paste0(getwd(), "/results/otter_trawl/"),
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
    return(NULL)
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
  model_settings$status[irow] <- ifelse(test = is.null(fit) == T , 
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
    
    ## Deviance and AIC
    # model_settings$deviance[irow] <- fit$Report$deviance
    model_settings$aic[irow] <- fit$parameter_estimates$AIC
    
  } else (model_settings[irow, c("status", "rrmse", "max_grad", "aic")] <- NA)
  
  rm(fit)
}

dyn.unload(paste0(getwd(), "/results/otter_trawl//", vast_cpp_version, ".dll"))

##################################################
####   Save
##################################################
write.csv(model_settings, 
          file = paste0("results/otter_trawl/model_settings.csv"),
          row.names = F)

lapply(X = split.data.frame(x = model_settings, 
                            f = model_settings$common_name),
       FUN = function(x) {
         temp_df <- subset(x, max_grad < 1e-4)
         min_idx <- which.min(temp_df$aic)
         
         return(temp_df[min_idx, ])
       })

model_settings <- read.csv(paste0("results/otter_trawl/model_settings.csv"))

##################################################
####   Loop over species and fit the model with "best" field configurations
##################################################

for (ispp in spp_list) { ## Loop over species -- start
  
  ## Subset data input for species ispp
  data_geostat_subset <- subset(x = data_geostat, 
                                subset = spp == ispp)
  
  ## Subset the model runs that converged with a low enough maximum gradient 
  sub_df <- subset(x = model_settings, 
                   subset = common_name == ispp & max_grad < 1e-4)
  
  ## If there is at least one converged model, run the model with the lowest AIC
  if (nrow(sub_df) > 0) {
    
    # best_model_idx <- which.max(1 - sub_df$deviance/max(sub_df$deviance))
    
    best_model_idx <- which.min(sub_df$aic)
    field_config <- unlist(sub_df[best_model_idx, 
                                  c("Omega1", "Epsilon1", 
                                    "Omega2", "Epsilon2")])
    settings$FieldConfig <- field_config
    
    ## Fit model
    result_dir <- paste0(getwd(), "/results/otter_trawl/", ispp)
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
    
    # dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
    # dyn.load("C:/Users/zack.oyafuso/Work/GitHub/Arctic_GF_OM/results/otter_trawl/Alaska plaice/VAST_v12_0_0.dll")
    ## Diagnostics
    diagnostics <- plot(fit, working_dir = paste0(result_dir, "/"))
    
    ## Save fit and diagnostics
    saveRDS(fit, paste0(result_dir, "/fit_full.rds")) # save all outputs locally
    fit <- fit[c("parameter_estimates", "data_frame", "data_list", "Report")] # partial output to sync to remote
    save(list = "fit", file = paste0(result_dir, "/fit.RData"))
    save(list = "diagnostics", file = paste0(result_dir, "/diagnostics.RData"))
    
    dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
    
  }
  
} ## Loop over species -- end
