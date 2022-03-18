#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       VAST Modelling, General Code for Chukchi
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   1990 and 2012 Chukchi otter trawl data
##                2012, 2017, 2019 Chukchi beam trawl data
##                No covariates
## 
## Notes:         Versions to use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Software settings ----
##    Check that versions of R and relevant packages are consistent
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
R_version <- "R version 4.0.2 (2020-06-22)"
VAST_cpp_version <- "VAST_v13_1_0"
pck_version <- c("VAST" = "3.9.0", "FishStatsUtils" = "2.11.0", 
                 "Matrix" = "1.4-0", "TMB" = "1.7.22", "DHARMa" = "0.4.5")

{
  if(sessionInfo()$R.version$version.string == R_version) 
    message(paste0(sessionInfo()$R.version$version.string, 
                   " IS CONSISTENT with the 2022 TOR."))
  
  if(!sessionInfo()$R.version$version.string == R_version) 
    message(paste0(sessionInfo()$R.version$version.string, 
                   " IS NOT CONSISTENT with the 2022 TOR. ",
                   "Please update R version to ", R_version))
  
  for (pck in 1:length(pck_version)) {
    temp_version <- packageVersion(pkg = names(pck_version)[pck])
    
    if(temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") IS CONSISTENT",
                     " with the 2022 TOR."))
    
    if(!temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") IS NOT CONSISTENT",
                     " with the 2022 TOR. Please update the '", 
                     names(pck_version)[pck], "' package to ", 
                     pck_version[pck]))
  }
  rm(pck, temp_version)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(VAST)
library(googledrive)
library(purrr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Authorize Google Drive ----
##   Create folder in google drive where the VAST output will be saved. 
##   This is done because VAST outputs and the simulated density objects are
##   often very large (i.e., > 100 MB) and go over github's data limits. This
##   makes it very bulky to push and pull these outputs over github.
##
##   For downstream analyses, VAST output is to be downloaded from google drive. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
googledrive::drive_deauth()
googledrive::drive_auth() 
1

if( nrow(googledrive::drive_get(path = "Oyafuso_Chukchi_VAST")) == 0 ) {
  google_folder <- googledrive::drive_mkdir("Oyafuso_Chukchi_VAST")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   VAST Model Settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

settings <- FishStatsUtils::make_settings(
  Version = VAST_cpp_version,
  n_x = 200,
  Region = "chukchi_sea", 
  purpose = "index2",
  ObsModel = c(2, 1),
  max_cells = Inf,
  use_anisotropy = TRUE, 
  Options = c('SD_site_logdensity' = FALSE, 'Calculate_Range' = TRUE,
              'Calculate_effective_area' = FALSE, 'Calculate_Cov_SE' = FALSE,
              'Calculate_Synchrony' = FALSE, 'Calculate_proportion' = FALSE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Main Loop ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Loop over gears: start ---- 
for (igear in c("beam", "otter")[2]) { 
  
  ## Create a folder in google drive for gear igear if it doesn't exist
  gear_folder <- paste0("Oyafuso_Chukchi_VAST", "/", igear, "/") 
  if(nrow(googledrive::drive_get(path = paste0("Oyafuso_Chukchi_VAST/", 
                                               igear, "/"))) == 0) {
    googledrive::drive_mkdir(path = "Oyafuso_Chukchi_VAST/", 
                             name = igear, 
                             overwrite = TRUE)
  }
  
  ## Species list
  spp_list <- 
    gsub(x = dir(path = paste0("data/fish_data/",
                               switch(igear,
                                      otter = "AK_BTS_OtterAndBeam",
                                      beam = "2017_2019_Beam"),
                               "/data_long_by_taxa/")), 
         pattern = ".csv",
         replacement = "")
  
  ## Loop over species: start ---- 
  for (ispp in spp_list[2] ) { 
    
    ## Create species folder locally
    result_dir <- paste0(getwd(), "/results/chukchi_", igear, 
                         "/vast_fits/", ispp, "/")
    if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    
    ## CPUE Data for species ispp
    AKBTS_data <- read.csv(paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                  "data_long_by_taxa/", 
                                  ispp, ".csv"))
    ierl_beam_data <- read.csv(paste0("data/fish_data/2017_2019_Beam/",
                                      "data_long_by_taxa/", 
                                      ispp, ".csv"))
    
    ## subset gear igear
    spp_data <- subset(x = AKBTS_data, subset = gear == igear,
                       select = c(year, lon, lat, cpue_kg_km2))
    
    if (igear == "beam") { #add ierl beam data
      spp_data <- rbind(spp_data,
                        subset(x = ierl_beam_data,
                               select = c(year, lon, lat, cpue_kg_km2)))  
    }
    
    ## Data input for VAST
    data_geostat <- data.frame(
      region = "chukchi",
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
    model_settings <- expand.grid( region = "chukchi",
                                   gear = igear,
                                   common_name = ispp,
                                   Omega1 = c(0, "IID"),
                                   Epsilon1 = c(0, "IID"),
                                   Omega2 = c(0, "IID"),
                                   Epsilon2 = c(0, "IID"),
                                   obs_model = c("PosLink", 
                                                 "Log_Delta", 
                                                 "Gamma_Delta"),
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
      
      ## Set observation model
      settings$ObsModel <- switch(model_settings$obs_model[irow], 
                                  "PosLink" = c(2, 4),
                                  "Log_Delta" = c(1, 0), 
                                  "Gamma_Delta" = c(2, 0))
      
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
      
      settings$ObsModel <- switch(sub_df$obs_model[best_model_idx], 
                                  "PosLink" = c(2, 4),
                                  "Log_Delta" = c(1, 0), 
                                  "Gamma_Delta" = c(2, 0))
      
      ## Fit model with best settings ----
      ## Turn on bias correction
      settings$bias.correct <- TRUE
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
      grid_df <- data.frame()
      for (itime in sort(unique(data_geostat$Year))) {
        grid_df <-
          rbind(grid_df,
                data.frame(spp = ispp,
                           Year = rep(itime, nrow(chukchi_sea_grid)),
                           Catch_KG = mean(data_geostat$Catch_KG),
                           AreaSwept_km2 = chukchi_sea_grid[, "Area_in_survey_km2"],
                           Lat = chukchi_sea_grid[, "Lat"],
                           Lon = chukchi_sea_grid[, "Lon"],
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
        "Parameters" = ParHat)
      
      ## Save local copy of the fit
      saveRDS(fit_sim, paste0(result_dir, 
                              "simulated_densities/fit_sim_full.rds"))
      
      ## Simulate Density ----
      ## Simulate_data() function produces simulated biomasses, and we 
      ## divide by the cell areas to get density.
      sim_data <- array(data = NA,
                        dim = c(nrow(chukchi_sea_grid),
                                length(unique(data_geostat$Year)),
                                1000))
      
      for (isim in 1:1000) {
        Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim,
                                              type = 1,
                                              random_seed = isim)
        
        sim_data[, ,isim] <- matrix(Sim1$b_i[pred_TF == 1] / 
                                      chukchi_sea_grid[, "Area_in_survey_km2"],
                                    nrow = nrow(chukchi_sea_grid))
        
        if(isim%%50 == 0) print(paste("Done with", ispp, "Iteration", isim))
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
      
      ## Upload to Google Drive ----
      
      ## Create species folder in google, Upload files related to the VAST fit
      if(nrow(googledrive::drive_get(
        path = paste0(gear_folder, ispp, "/"))) == 0) {
        spp_folder <- googledrive::drive_mkdir(path = gear_folder, 
                                               name = ispp)
      }
      
      googledrive::with_drive_quiet(
        files <- purrr::map(paste0(result_dir, 
                                   c("Kmeans_knots-200.RData", 
                                     "model_settings.csv", 
                                     "packageDescription.txt", 
                                     "parameter_estimates.RData", 
                                     "parameter_estimates.txt", 
                                     "settings.txt",
                                     "fit.RData", 
                                     "fit_full.rds")), 
                            ~ drive_upload(.x, path = spp_folder)))
      
      ## Create diagnostics folder in Google Drive, Upload Diagnostics
      if(nrow(googledrive::drive_get(
        path = paste0(gear_folder, ispp, "/diagnostics/"))) == 0) {
        diagnostics_folder <- 
          googledrive::drive_mkdir(paste0("Oyafuso_Chukchi_VAST/",
                                          igear, "/", ispp, 
                                          "/diagnostics"))
      }
      
      googledrive::with_drive_quiet(
        files <- purrr::map(dir(paste0(result_dir, "diagnostics/"), 
                                full.names = T), 
                            ~ drive_upload(.x, path = diagnostics_folder )))
      
      ## Create simulated_densities folder in Google Drive, upload 
      ## simulated densities
      if(nrow(googledrive::drive_get(
        path = paste0(gear_folder, ispp, "/simulated_densities/"))) == 0) {
        simualtion_folder <- 
          googledrive::drive_mkdir(paste0("Oyafuso_Chukchi_VAST/",
                                          igear, "/", ispp, 
                                          "/simulated_densities"))
      }
      
      googledrive::with_drive_quiet(
        files <- purrr::map(dir(paste0(result_dir, "simulated_densities/"),
                                full.names = T), 
                            ~ drive_upload(.x, path = simualtion_folder )))
      
    }  
  } ## Loop over species: end
} ## Loop over gears: end 

