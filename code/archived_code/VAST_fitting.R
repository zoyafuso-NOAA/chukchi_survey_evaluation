###############################################################################
## Project:       VAST Modelling, General Code for Chukchi
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   1990 and 2012 Chukchi otter trawl data
##                2012, 2017, 2019 Chukchi beam trawl data
##                No covariates
## 
## Notes:         Versions to use
rm(list = ls())
r_version <- "R version 4.0.2 (2020-06-22)"
vast_version <- "3.6.1"
fishstatsutils_version <- "2.8.0"
vast_cpp_version <- "VAST_v12_0_0"
###############################################################################

##################################################
####  Check that versions of R and relevant packages are consistent 
##################################################
ifelse(test = sessionInfo()$R.version$version.string == r_version &
         packageVersion("VAST") == vast_version &
         packageVersion("FishStatsUtils") == fishstatsutils_version &
         FishStatsUtils::get_latest_version() == vast_cpp_version,
       yes = "versions are good to go", 
       no = "check your versions")

##################################################
####  Import Libraries
##################################################
# library(devtools)
# devtools::install_local("C:/Users/zack.oyafuso/Downloads/FishStatsUtils-2.8.0")
# devtools::install_github("James-Thorson-NOAA/VAST@3.6.1") #for a specific VAST package version
library(VAST)

##################################################
####  Import interpolation grid
##################################################
chukchi_grid <- read.csv(file = paste0("data/spatial_data/",
                                       "BS_Chukchi_extrapolation_grids/",
                                       "ChukchiThorsonGrid.csv"))
chukchi_grid$Area_km2 <- chukchi_grid$Shape_Area / 1000 / 1000

##################################################
####  Synthesize a main cpue dataframe
##################################################
aksc_data <- read.csv(file = paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                    "AK_BTS_Arctic_processed_long.csv"))
aksc_data <- subset(x = aksc_data, 
                    subset = common_name %in% c("yellowfin sole", "snow crab",
                                                "saffron cod",
                                                "Bering flounder",
                                                "Arctic cod", "Alaska plaice"))
chukchi_beam_data <- read.csv(file = paste0("data/fish_data/2017_2019_Beam/",
                                            "ierl_data_processed.csv"))
chukchi_beam_data <- subset(x = chukchi_beam_data,
                            subset = common_name %in% c("yellowfin sole",
                                                        "snow crab",
                                                        "saffron cod",
                                                        "Bering flounder",
                                                        "Arctic cod"))

data_geostat <- data.frame( 
  gear = c(aksc_data$gear, chukchi_beam_data$gear),
  spp = as.factor(c(aksc_data$common_name, chukchi_beam_data$common_name)),
  Year = c(aksc_data$year, chukchi_beam_data$year),
  Catch_KG = c(aksc_data$cpue_kg_km2, chukchi_beam_data$cpue_kg_km2),
  AreaSwept_km2 = 1,
  Lat = c(aksc_data$lat, chukchi_beam_data$lat),
  Lon = c(aksc_data$lon, chukchi_beam_data$lon),
  stringsAsFactors = T)

##################################################
####   VAST Model Settings
##################################################
settings <- FishStatsUtils::make_settings( 
  Version = vast_cpp_version,
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
spp_list <- sort(unique(data_long$common_name))
model_settings <- 
  rbind(expand.grid(gear = "beam",
                    common_name = sort(unique(chukchi_beam_data$common_name)),
                    Omega1 = 0:1,
                    Omega2 = 0:1,
                    Epsilon1 = 0:1,
                    Epsilon2 = 0:1,
                    stringsAsFactors = FALSE),
        expand.grid(gear = "otter",
                    common_name = sort(unique(aksc_data$common_name)),
                    Omega1 = 0:1,
                    Omega2 = 0:1,
                    Epsilon1 = 0:1,
                    Epsilon2 = 0:1,
                    stringsAsFactors = FALSE))

model_settings[, c("status", "max_grad", "rrmse", "aic")] <- NA

##################################################
####   Create result directory if not created already
##################################################
if(!dir.exists("results/otter_trawl/vast_fits/")) 
  dir.create("results/otter_trawl/vast_fits/")
if(!dir.exists("results/beam_trawl_2012_2019/vast_fits/"))
  dir.create("results/beam_trawl_2012_2019/vast_fits/")

##################################################
####   Model fit
##################################################
for (irow in 1:nrow(model_settings)) {
  
  result_dir <- paste0(getwd(), "/results/", 
                       ifelse(test = model_settings$gear[irow] == "otter",
                              yes = "otter_trawl/vast_fits/",
                              no = "beam_trawl_2012_2019/vast_fits/"))
  
  data_geostat_subset <- 
    subset(x = data_geostat, 
           subset = spp == model_settings$common_name[irow] &
             gear == model_settings$gear[irow])
  
  settings$FieldConfig <- unlist(model_settings[irow, c("Omega1", "Epsilon1",
                                                        "Omega2", "Epsilon2")])
  
  fit <- tryCatch( {FishStatsUtils::fit_model( 
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
    "input_grid" = chukchi_grid
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
  
  dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
}



##################################################
####   Save
##################################################
for (igear in c("beam", "otter")) {
  write.csv(x = subset(x = model_settings, subset = gear == igear), 
            file = paste0("results/",
                          ifelse(test = igear == "otter",
                                 yes = "otter_trawl/vast_fits/",
                                 no = "beam_trawl_2012_2019/vast_fits/"),
                          "model_settings.csv"),
            row.names = F)
  
}

##################################################
####   Loop over species and fit the model with "best" field configurations
##################################################
for (igear in c("beam", "otter")) {  ## Loop over gear
  
  spp_list <- paste(unique(data_geostat$spp[data_geostat$gear == igear]))
  n_years <- length(unique(data_geostat$Year[data_geostat$gear == igear]))
  years <- sort(unique(data_geostat$Year[data_geostat$gear == igear]))
  
  for (ispp in spp_list) { ## Loop over species -- start
    
    ## Subset data input for species ispp
    data_geostat_subset <- subset(x = data_geostat, 
                                  subset = gear == igear & spp == ispp)
    
    ## Subset the model runs that converged with a low enough maximum gradient 
    sub_df <- subset(x = model_settings, 
                     subset = common_name == ispp 
                     & gear == igear
                     & max_grad < 1e-4)
    
    ## If there is at least one converged model, run the model with the lowest AIC
    if (nrow(sub_df) > 0) {
      
      best_model_idx <- which.min(sub_df$aic)
      field_config <- unlist(sub_df[best_model_idx, 
                                    c("Omega1", "Epsilon1", 
                                      "Omega2", "Epsilon2")])
      settings$FieldConfig <- field_config
      
      ## Fit model directory, copy VAST cpp files
      result_dir <- paste0(getwd(), "/results/", 
                           ifelse(test = igear == "otter",
                                  yes = "otter_trawl/vast_fits/",
                                  no = "beam_trawl_2012_2019/vast_fits/"),
                           ispp, "/")
      if(!dir.exists(result_dir)) dir.create(result_dir)
      
      ##################################################
      ####   Fit the model and save output
      ##################################################
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
        "input_grid" = chukchi_grid)
      
      ## Diagnostics
      diagnostics <- plot(fit, working_dir = paste0(result_dir, "/"))
      save(list = "diagnostics", 
           file = paste0(result_dir, "/diagnostics.RData"))
      
      ## Save fit and diagnostics
      saveRDS(fit, paste0(result_dir, "/fit_full.rds"))
      
      ## partial output to sync to remote
      fit <- fit[c("parameter_estimates", "data_frame", "data_list", "Report")] 
      save(list = "fit", file = paste0(result_dir, "/fit.RData"))
      
      dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
      
      ##################################################
      ####   Refit the model with grid locations attached to data in order
      ####       to simulate density values
      ##################################################
      ## Prediction Grid: df of the grid to simulate data onto
      grid_df <- data.frame()
      for (itime in years) {
        grid_df <- 
          rbind(grid_df,
                data.frame(spp = ispp,
                           Year = rep(itime, nrow(chukchi_grid)),
                           Catch_KG = mean(data_geostat_subset$Catch_KG),
                           AreaSwept_km2 = chukchi_grid$Area_km2,
                           Lat = chukchi_grid$Lat,
                           Lon = chukchi_grid$Lon,
                           stringsAsFactors = T)
          )
      }
      
      ###################################################
      ## Add New Points: set catch to NAs?
      ###################################################
      data_geostat_with_grid <- rbind(data_geostat_subset[, names(grid_df)],
                                      grid_df)
      
      pred_TF <- rep(1, nrow(data_geostat_with_grid))
      pred_TF[1:nrow(data_geostat_subset)] <- 0
      
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
        "input_grid" = chukchi_grid, 
        "PredTF_i" = pred_TF)
      
      ## Save local copy of the fit
      saveRDS(fit_sim, paste0(result_dir, "/fit_sim_full.rds"))
      
      ##################################################
      ####   Simulate 500 iterations of densities
      ####   Simulat_data() function produces simulated biomasses, and we 
      ####      divide by the cell areas to calculate simulated desniteis
      ##################################################
      sim_data <- array(data = NA, 
                        dim = c(nrow(chukchi_grid), n_years, 500))
      
      for (isim in 1:500) {
        Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim, 
                                              type = 1, 
                                              random_seed = isim)
        sim_data[, ,isim] <- matrix(Sim1$b_i[pred_TF == 1] / chukchi_grid$Area_km2, 
                                    nrow = nrow(chukchi_grid))
        if(isim%%50 == 0) print(paste("Done with", ispp, "Iteration", isim))
      }
      
      ##################################################
      ####   Save simulated densities by year so they can be pushed to github
      ##################################################
      for (iyear in 1:n_years) {
        assign(value = sim_data[, iyear, ], 
               x = paste0("sim_data_", years[iyear]) )
        save(list = paste0("sim_data_", years[iyear]), 
             file = paste0(result_dir, "/sim_data_", 
                           years[iyear], ".RData") )
      }
      
      ## partial output to sync to remote
      fit_sim <- fit_sim[c("parameter_estimates", "data_frame", 
                           "data_list", "Report")] 
      save(list = "fit_sim", 
           file = paste0(result_dir, "/fit_sim.RData"))
      
      dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
      
    }
  }  ## Loop over species -- end
}


