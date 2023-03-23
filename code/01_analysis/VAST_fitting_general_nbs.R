#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Spatiotemporal distribution modelling using VAST
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Lewis Barnett (lewis.barnett@noaa.gov)
## Description:   1982 to 2022 NBS otter trawl data analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(VAST)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   VAST Model Settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
settings <- FishStatsUtils::make_settings(
  n_x = 200,
  Region = "northern_bering_sea", 
  purpose = "index2",
  ObsModel = c(2, 1),
  max_cells = Inf,
  use_anisotropy = FALSE, 
  Options = c('SD_site_logdensity' = FALSE, 'Calculate_Range' = FALSE,
              'Calculate_effective_area' = FALSE, 'Calculate_Cov_SE' = FALSE,
              'Calculate_Synchrony' = FALSE, 'Calculate_proportion' = FALSE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Main Loop ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Species list 
spp_list <- 
  gsub(x = dir(path = paste0("data/fish_data/AK_BTS_OtterAndBeam/data_long_by_taxa_nbs/")), 
       pattern = ".csv",
       replacement = "")

## Loop over species: start ---- 
for (ispp in spp_list) { 
  
  ## Create species folder locally
  result_dir <- paste0(getwd(), "/results/nbs_otter/vast_fits/", ispp, "/")
  if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = TRUE)
  
  ## CPUE Data for species ispp
  AKBTS_data <- read.csv(file = paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                       "data_long_by_taxa_nbs/", ispp, ".csv"))
  
  spp_data <- select(AKBTS_data, year, lon, lat, cpue_kg_km2)
  
  ## Data input for VAST
  data_geostat <- data.frame(
    region = "northern_bering_sea",
    spp = ispp,
    Year = spp_data$year,
    Catch_KG = spp_data$cpue_kg_km2,
    AreaSwept_km2 = 1,
    Lat = spp_data$lat,
    Lon = spp_data$lon,
    stringsAsFactors = T)
  
  ## Model settings ----
  ## We run VAST under varying configurations of the Spatial Field (Omega1
  ## and Omega2), Spatiotemporal Field (Epsilon1 and Epsilon2), and 
  ## Observation Model (obs_model). Each scenario's settings are kept in 
  ## the model_settings df.
  model_settings <- expand.grid(region = "northern_bering_sea",
                                common_name = ispp,
                                Omega1 = c("IID"),
                                Epsilon1 = c(0, "IID"),
                                Omega2 = c("IID"),
                                Epsilon2 = c(0, "IID"),
                                stringsAsFactors = FALSE)
  model_settings[, c("status", "max_grad", "rrmse", "aic")] <- NA
  settings$bias.correct <- FALSE ## Turn off during initial runs
  
  ## Loop over Model Configs: start ---- 
  for (irow in 1:nrow(model_settings)) { 
    
    ## Set Spatial/Spatiotemporal fields
    settings$FieldConfig <- unlist(model_settings[irow, c("Omega1", 
                                                          "Epsilon1", 
                                                          "Omega2", 
                                                          "Epsilon2")])
    
    ## Fit model. The tryCatch function makes sure that if the model does
    ## not converge, it does not break the loop.
    fit <- tryCatch( {FishStatsUtils::fit_model( 
      "settings" = settings,
      "working_dir" = result_dir,
      "Lat_i" = data_geostat[, "Lat"],
      "Lon_i" = data_geostat[, "Lon"],
      "t_i" = data_geostat[, "Year"],
      "b_i" = data_geostat[, "Catch_KG"],
      "a_i" = data_geostat[, "AreaSwept_km2"],
      "getJointPrecision" = TRUE,
      "newtonsteps" = 1,
      "test_fit" = F)
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
      message( paste0("\nProcessed model ", irow, 
                      " of ", nrow(model_settings), " for ", ispp, ": ", 
                      model_settings$gear[irow], ", ",
                      model_settings$common_name[irow],
                      ", with settings\n",
                      model_settings$obs_model[irow],
                      ", Omega1 = ", model_settings$Omega1[irow],
                      " Omega2 = ", model_settings$Omega2[irow],
                      " Epsilon1 = ", model_settings$Epsilon1[irow],
                      " Epsilon2 = ", model_settings$Epsilon2[irow]) )
    })    
    
    ## Check whether the model converges
    model_settings$status[irow] <- 
      ifelse(test = is.null(fit) == T | 
               is.null(fit$parameter_estimates$max_gradient) , 
             yes = "no_convergence",
             no = "confirm_gradient")
    
    ## If the model runs ----
    ## Collect the max gradient, rrmse of predictions, and aic
    if (model_settings$status[irow]  == "confirm_gradient") {
      
      ## Maximum Gradient
      model_settings$max_grad[irow] <- fit$parameter_estimates$max_gradient
      
      ## RRMSE of density predictions
      obs_cpue <- with(data_geostat,  Catch_KG / AreaSwept_km2 )
      pred_cpue <- fit$Report$D_i
      rrmse <- sqrt(mean((obs_cpue - pred_cpue)^2)) / mean(obs_cpue)
      
      model_settings$rrmse[irow] <- round(rrmse, 3)
      
      ## Deviance and AIC
      model_settings$aic[irow] <- fit$parameter_estimates$AIC
      
    } else (model_settings[irow, c("status", "rrmse", 
                                   "max_grad", "aic")] <- NA)
    
    rm(fit)
    
  } ## Loop over Model Configs: end
  
  ## Save model settings
  write.csv(x = model_settings, 
            file = paste0(result_dir, "model_settings.csv"),
            row.names = F)
  
  ## Find model with best settings ----
  ## Subset the model runs that converged with a low enough maximum gradient 
  sub_df <- subset(x = model_settings, subset = max_grad < 1e-4)
  
  if (nrow(sub_df) > 0) {
    best_model_idx <- which.min(sub_df$aic)
    field_config <- unlist(sub_df[best_model_idx, 
                                  c("Omega1", "Epsilon1", 
                                    "Omega2", "Epsilon2")])
    settings$FieldConfig <- field_config
    
    ## Fit model with best settings ----
    ## Turn on bias correction
    settings$bias.correct <- FALSE
    fit <- FishStatsUtils::fit_model( 
      "settings" = settings,
      "working_dir" = result_dir,
      "Lat_i" = data_geostat[, "Lat"],
      "Lon_i" = data_geostat[, "Lon"],
      "t_i" = data_geostat[, "Year"],
      "b_i" = data_geostat[, "Catch_KG"],
      "a_i" = data_geostat[, "AreaSwept_km2"],
      "getJointPrecision" = TRUE,
      "newtonsteps" = 1,
      "test_fit" = FALSE)
    
    ## Diagnostics ----
    if(!dir.exists(paste0(result_dir, "diagnostics/")))
      dir.create(paste0(result_dir, "diagnostics/"))
    
    diagnostics <- plot(fit, working_dir = paste0(result_dir, "diagnostics/"))
    save(list = "diagnostics", 
         file = paste0(result_dir, "diagnostics/diagnostics.RData"))
    
    ## Save fit and diagnostics
    saveRDS(fit, paste0(result_dir, "fit_full.rds"))
    
    ## Save parameter estimates
    ParHat <- fit$ParHat
    
    ## partial output to sync to remote
    fit <- fit[c("parameter_estimates", "data_frame", "data_list", "Report")] 
    save(list = "fit", file = paste0(result_dir, "fit.RData"))
    
    ## Refit VAST w/ grid locs ----
    ## This is a workaround to simulate density values across each grid cell
    if(!dir.exists(paste0(result_dir, "simulated_densities/")))
      dir.create(paste0(result_dir, "simulated_densities/"))
    
    ## Prediction Grid: first, create a dataframe of locations, grid cell
    ## areas and dummy catches (set to the mean) for each year.
    nbs_grid <- read.csv("data/spatial_data/BS_Chukchi_extrapolation_grids/nbs_2022.csv")
    grid_df <- data.frame()
    for (itime in sort(unique(data_geostat$Year))) {
      grid_df <-
        rbind(grid_df,
              data.frame(spp = ispp,
                         Year = rep(itime, nrow(nbs_grid)),
                         Catch_KG = mean(data_geostat$Catch_KG),
                         AreaSwept_km2 = nbs_grid[, "Area_KM2"],
                         Lat = nbs_grid[, "lat"],
                         Lon = nbs_grid[, "lon"],
                         stringsAsFactors = T)
        )
    }
    
    ## Second, append grid locations to data_geostat
    data_geostat_with_grid <- rbind(data_geostat[, names(grid_df)],
                                    grid_df)
    
    ## Third, specifiy pred_TF, a binary vector with values:
    ## 0: use in model to estimate VAST parameters
    ## 1: ignore in model fitting but still predict densities onto
    ## This is our workaround to be able to simulate densities for each
    ## grid cell.
    pred_TF <- rep(1, nrow(data_geostat_with_grid))
    pred_TF[1:nrow(data_geostat)] <- 0
    
    ## Fourth: refit model with expanded grid data. Use parameter estimate
    ## from the original fit (ParHat) as a shortcut.
    fit_sim <- FishStatsUtils::fit_model(
      "settings" = settings,
      "working_dir" = result_dir,
      "Lat_i" = data_geostat_with_grid[, "Lat"],
      "Lon_i" = data_geostat_with_grid[, "Lon"],
      "t_i" = data_geostat_with_grid[, "Year"],
      "b_i" = data_geostat_with_grid[, "Catch_KG"],
      "a_i" = data_geostat_with_grid[, "AreaSwept_km2"],
      "getJointPrecision" = TRUE,
      "newtonsteps" = 1,
      "test_fit" = F,
      
      "PredTF_i" = pred_TF,
      "Parameters" = ParHat
    )
    
    ## Save local copy of the fit
    saveRDS(fit_sim, paste0(result_dir, 
                            "simulated_densities/fit_sim_full.rds"))
    
    ## Simulate Density ----
    ## Simulate_data() function produces simulated biomasses, and we 
    ## divide by the cell areas to get density.
    sim_data <- array(data = NA,
                      dim = c(nrow(nbs_grid),
                              length(unique(data_geostat$Year)),
                              1000))
    
    for (isim in 1:1000) {
      Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim,
                                            type = 1,
                                            random_seed = isim)
      
      sim_data[, ,isim] <- matrix(Sim1$b_i[pred_TF == 1] / 
                                    nbs_grid[, "Area_KM2"],
                                  nrow = nrow(nbs_grid))
      
      if(isim%%10 == 0) print(paste("Done with", ispp, "Iteration", isim))
    }
    
    ## Save simulated densities by year locally
    for (iyear in 1:length(unique(data_geostat$Year))) {
      
      obj_name <- paste0("sim_data_",
                         sort(unique(data_geostat$Year))[iyear],
                         "_simtype1")
      
      assign(value = sim_data[, iyear, ],
             x = obj_name )
      save(list = obj_name,
           file = paste0(result_dir, "simulated_densities/", 
                         obj_name, ".RData") )
    }
    
    ## Save partial output 
    fit_sim <- fit_sim[c("parameter_estimates", "data_frame",
                         "data_list", "Report")]
    save(list = "fit_sim",
         file = paste0(result_dir, "simulated_densities/fit_sim.RData"))
    
  }  
} ## Loop over species: end


