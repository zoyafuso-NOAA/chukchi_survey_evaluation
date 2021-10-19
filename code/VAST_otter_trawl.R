###############################################################################
## Project:       VAST Modelling, Chukchi Sea
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   1990 and 2012 Chukchi data
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries
##################################################
library(tidyr)
library(reshape)

library(raster)
library(sp)
library(rgdal)

# library(devtools)
# devtools::install_github("James-Thorson-NOAA/VAST@3.7.1") #for a specific VAST package version
library(VAST)

##################################################
####  Import Spatial data
##################################################
# AK <- rgdal::readOGR(paste0(wd, "spatial_data/land_shapefiles/AKland.shp"))
# AK <- sp::spTransform(AK, CRSobj = crs("+proj=longlat") )

##################################################
####  Import Data
##################################################
<<<<<<< Updated upstream
rb_data <- read.csv("data/fish_data/otter_trawl/AK_BTS_Arctic_processed_long.csv")
ierl_data <- read.csv("data/fish_data/2017_2019_Beam/ierl_data_processed.csv")
=======
rb_data <- read.csv("data/fish_data/otter_trawl/AK_BTS_Arctic_processed.csv")
# ierl_data <- read.csv("data/fish_data/2017_2019_Beam/ierl_data_processed.csv")
>>>>>>> Stashed changes

extrapolation_grid <- read.csv(paste0("data/spatial_data/",
                                      "BS_Chukchi_extrapolation_grids/",
                                      "ChukchiThorsonGrid.csv"))
extrapolation_grid$Area_km2 <- extrapolation_grid$Shape_Area / 1000 / 1000

<<<<<<< Updated upstream
data_long <- rbind(rb_data[, c("year", "area_swept_km2", "lat", "lon", 
                               "common_name", "catch_kg", "gear")],
                   ierl_data[, c("year", "area_swept_km2", "lat", "lon", 
                                 "common_name", "catch_kg", "gear")])
=======
data_long <- rb_data #rbind(rb_data, ierl_data[, names(rb_data)])
>>>>>>> Stashed changes
data_long$cpue <- data_long$catch_kg / data_long$area_swept_km2

spp_list <- c("Alaska plaice", "Arctic cod", "Bering flounder", "saffron cod",
              "walleye pollock")

data_long <- subset(x = data_long,
                    subset = gear == "otter" & 
                      common_name %in% spp_list &
                      year %in% c(1990, 2012))

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
VAST_Version <- "VAST_v13_1_0"
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
              'Calculate_Synchrony' = FALSE, 'Calculate_proportion' = FALSE,
              "report_additional_variables"=TRUE ))


##################################################
####   Model result objects
##################################################
<<<<<<< Updated upstream
model_settings <- expand.grid(common_name = spp_list,
                              Omega1 = 0:1,
                              Omega2 = 0:1,
                              Epsilon1 = 0:1,
                              Epsilon2 = 0:1,
=======
model_settings <- expand.grid(Omega1 = 0:1,
                              Omega2 = 0:1,
                              Epsilon1 = 0:1,
                              Epsilon2 = 0:1,
                              dens_covar = c("TRUE", "FALSE"),
                              common_name = spp_list,
>>>>>>> Stashed changes
                              stringsAsFactors = FALSE)
model_settings[, c("Beta1","Beta2")] <- 0

null_model_idx <- apply(X = model_settings[, c("Omega1", "Epsilon1",
                                               "Omega2", "Epsilon2")],
                        MARGIN = 1, 
                        FUN = function(x) all(x == 0))

model_settings[null_model_idx, c("Beta1","Beta2")] <- 3

model_settings[, c("status", "max_grad", "rrmse", "deviance")] <- NA

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
  settings$RhoConfig[c("Beta1","Beta2")] <- 
    unlist(model_settings[irow, c("Beta1","Beta2")])
  
  data_geostat_subset <- subset(x = data_geostat, 
                                subset = spp == model_settings$common_name[irow])
  
  fit <- tryCatch(
    {switch(paste0(model_settings$dens_covar[irow]),
            "TRUE" = FishStatsUtils::fit_model( 
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
              "input_grid" = extrapolation_grid,
              
              "X1_formula" =  "Catch_KG ~ Depth + Temp",
              "X2_formula" =  "Catch_KG ~ Depth + Temp",
              "covariate_data" = cbind(data_geostat_subset[, c("Lat",
                                                        "Lon",
                                                        "Catch_KG",
                                                        "Depth",
                                                        "Temp")],
                                       Year = NA)),
            
            "FALSE" = FishStatsUtils::fit_model( 
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
    )
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
  model_settings$status[irow] <- ifelse(test = (length(fit$Report) == 1) | (is.null(fit) == T) , 
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
    model_settings$deviance[irow] <- fit$Report$deviance
    model_settings$AIC[irow] <- fit$parameter_estimates$opt$AIC
    
  } else (model_settings[irow, c("status", "rrmse", "max_grad", "deviance")] <- NA)
}

## Unlink Dynamic library
dyn.unload(paste0(getwd(), "/results/otter_trawl//", VAST_Version, ".dll"))

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
         min_idx <- which.min(temp_df$rrmse)
         
         return(temp_df[min_idx, ])
       })

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
  
  ## If there is at least one converged model, run the model with the 
  ## highest deviance explained
  if (nrow(sub_df) > 0) {
    
    best_model_idx <- which.max(1 - sub_df$deviance/max(sub_df$deviance))
    
    field_config <- unlist(sub_df[best_model_idx, 
                                  c("Omega1", "Epsilon1", 
                                    "Omega2", "Epsilon2")])
    rho_config <-  unlist(sub_df[best_model_idx, 
                                 c("Beta1", "Beta2")])
    settings$FieldConfig <- field_config
    settings$RhoConfig[c("Beta1", "Beta2")] <- rho_config
    
    depth_in_model <- paste0(sub_df$dens_covar[which.min(sub_df$rrmse)])
    
    ## Fit model
    result_dir <- paste0(getwd(), "/results/otter_trawl/", ispp, "/")
    fit <- switch(depth_in_model,
                  "TRUE" = FishStatsUtils::fit_model( 
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
                    "input_grid" = extrapolation_grid,
                    
                    "X1_formula" =  "Catch_KG ~ Depth + I(Depth^2) + Temp + I(Temp^2)",
                    "X2_formula" =  "Catch_KG ~ Depth + I(Depth^2) + Temp + I(Temp^2)",
                    "covariate_data" = cbind(data_geostat_subset[, c("Lat",
                                                              "Lon",
                                                              "Catch_KG",
                                                              "Depth",
                                                              "Temp")],
                                             Year = NA)),
                  
                  "FALSE" = FishStatsUtils::fit_model( 
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
    )
    
    
    ## Diagnostics
    diagnostics <- plot(fit, working_dir = result_dir)
    
    ## Save fit and diagnostics
    fit <- fit[c("parameter_estimates", "data_frame", "data_list", "Report")]
    save(list = "fit", file = paste0(result_dir, "/fit.RData"))
    save(list = "diagnostics", file = paste0(result_dir, "/diagnostics.RData"))

  }
  
} ## Loop over species -- end

