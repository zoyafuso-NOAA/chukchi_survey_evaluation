###############################################################################
## Project:       Simulate Surveys Function
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################
## Import Libraries
##################################
library(sp)
library(rgeos)
library(raster)
library(rgdal)
library(viridis )

##################################
## Load input data and single-species survey outputs
## Load Chukchi grid
##################################
load("data/survey_opt_data/optimization_data.RData")
chukchi_grid <- read.csv(paste0("data/spatial_data/",
                                "BS_Chukchi_extrapolation_grids/",
                                "ChukchiThorsonGrid.csv"))
chukchi_grid$Area_km2 <- chukchi_grid$Shape_Area / 1000 / 1000
cell_area <- chukchi_grid$Area_km2
total_area <- sum(cell_area)

chukchi_mask <- 
  rgdal::readOGR("data/spatial_data/survey_boundary/CHUKCHI_2012.shp")
crs(chukchi_mask) <- aea_crs

##################################
## Calculate locations of systematic designs with different sampling efforts
##################################
{
  grid_pts_aea <- sp::SpatialPointsDataFrame(coords = grid_pts@coords,
                                             data = data.frame(id = 1:n_cells), 
                                             proj4string = aea_crs)
  temp_n <- 0
  temp_res <- 70000
  sys_settings <- NULL
  sys_idx <- list()
  
  png(filename = "presentations/results_11_30_2021/sys_survey_grids.png", 
      width = 4, height = 6, units = "in", res = 500)
  par(mfrow = c(3, 2), mar = c(0, 0, 0, 0))
  while (temp_n < 100) {
    
    ## Make grid based on latlon coords
    temp_grid <- sp::makegrid(x = grid_pts_aea, cellsize = temp_res)
    temp_grid <- sp::SpatialPoints(coords = coordinates(temp_grid), 
                                   proj4string = aea_crs)
    temp_grid <- crop(x = temp_grid, y = chukchi_mask)
    temp_n <- length(temp_grid)
    
    # plot(chukchi_mask); points(temp_grid) 
    
    grid_idx <- apply(X = gDistance(spgeom1 = grid_pts_aea, 
                                    spgeom2 = temp_grid, 
                                    byid = TRUE),
                      MARGIN = 1, 
                      FUN = which.min)
    
    plot(chukchi_mask, asp = 1)
    box()
    legend("bottomright", legend = paste0("n = ", temp_n), bty = "n", cex = 2)
    points(temp_grid,  pch = 16, cex = 1)
    
    sys_settings <- rbind(sys_settings, 
                          data.frame(n = temp_n,
                                     res = temp_res))
    
    temp_res <- temp_res - 5000
    sys_idx <- c(sys_idx, list(grid_idx))
  }
  
  dev.off()
}

##################################
##
##################################
igear = "otter"
true_index <- get(paste0("true_index_", igear))
true_index <- true_index[nrow(true_index), ]

# load(paste0("results/",
#             ifelse(test = igear == "otter",
#                    yes = "otter_trawl",
#                    no = "beam_trawl_2012_2019"),
#             "/survey_opt/SS/ss_knit_results_", igear, ".RData"))


# load(paste0("results/",
#             ifelse(test = igear == "otter",
#                    yes = "otter_trawl",
#                    no = "beam_trawl_2012_2019"),
#             "/survey_opt/MS/Run_",
#             ifelse(test = igear == "otter",
#                    yes = 14,
#                    no = 11),
#             "/result_list.RData"))

##################################
## Synthesize simulated densities
##################################
n_spp <-  get(paste0("n_spp_", igear))
spp_list <-  get(paste0("spp_list_", igear))
years <- switch(igear, 
                "beam" = c(2019),
                "otter" = c(2012))
n_years <- length(years)

ms_dens <- array(dim = c(n_spp, n_cells, 500), 
                 dimnames = list(spp_list, NULL, NULL))

for (ispp in 1:n_spp) {
  for (iyear in years){
    ## Load simulated densities
    load(paste0("results/",
                ifelse(test = igear == "otter",
                       yes = "otter_trawl",
                       no = "beam_trawl_2012_2019"),
                "/vast_fits/", spp_list[ispp], "/sim_data_", iyear, ".RData"))
    
    ms_dens[ispp, , ] <- get(paste0("sim_data_", iyear))
  }
}

##################################
## Simulate Simple Random and Systematic Designs at varying sampling efforts
##################################
srs_n <- seq(from = 40, to = 100, by = 10)
index_srs <- cv_srs <- rb_srs <- 
  array(dim = c(length(srs_n), n_spp, 500),
        dimnames = list(paste0("srs_n = ", srs_n), spp_list, NULL))

index_sys <- cv_sys <- rb_sys <- 
  array(dim = c(nrow(sys_settings), n_spp, 500),
        dimnames = list(paste0("sys_n = ", sys_settings$n), spp_list, NULL))

for (iter in 1:500) {
  
  ## Simple Random Design
  for (isample in 1:length(srs_n)) {
    
    set.seed(iter * 234)
    srs_sample_idx <- sample(x = 1:n_cells, 
                             size = srs_n[isample])
    srs_sample <- ms_dens[, srs_sample_idx, iter]
    srs_tau <- total_area * rowMeans(srs_sample)
    
    srs_sd <- sqrt(apply(X = srs_sample, MARGIN = 1, FUN = var) * 
                     total_area^2 / srs_n[isample])
    
    index_srs[isample, , iter] <- srs_tau / 1000
    rb_srs[isample, , iter] <- 100 * (srs_tau / 1000 - true_index) / true_index
    cv_srs[isample, , iter] <-  srs_sd / srs_tau
  }
  
  ## Systematic Design
  for (isample in 1:nrow(sys_settings)) {
    
    Y <- ms_dens[, sys_idx[[isample]], iter]
    # n_sys <- sum(cell_area[sys_idx[[isample]]])
    sys_var <- apply(X = Y, MARGIN = 1, FUN = var)
    sys_tau <- rowMeans(Y) * total_area
    
    sys_sd <- sqrt(total_area^2 * sys_var / sys_settings$n)
    cv_sys[isample, , iter] <- sys_sd / sys_tau
    
    index_sys[isample, , iter] <- sys_tau / 1000
    rb_sys[isample, , iter] <- 100 * (sys_tau / 1000 - true_index) / true_index
  }
}

## Plots
true_cv_srs <- sweep(x = apply(X = index_srs, MARGIN = 1:2, FUN = sd), 
                     MARGIN = 2, 
                     STATS = true_index, 
                     FUN = "/")

rrmse_cv_srs <- sqrt(apply(X = sweep(x = cv_srs, 
                                     MARGIN = 1:2, 
                                     STATS = true_cv_srs, 
                                     FUN = "-") ^ 2, ## Square Error
                           MARGIN = 1:2, 
                           FUN = mean) ## Mean Square Error
) /## Root Mean Square Error
  apply(X = cv_srs, MARGIN = 1:2, FUN = mean) ## Normalize by mean(cvs)


true_cv_sys <- sweep(x = apply(X = index_sys, MARGIN = 1:2, FUN = sd), 
                     MARGIN = 2, 
                     STATS = true_index, 
                     FUN = "/")

rrmse_cv_sys <- sqrt(apply(X = sweep(x = cv_sys, 
                                     MARGIN = 1:2, 
                                     STATS = true_cv_sys, 
                                     FUN = "-") ^ 2, ## Square Error
                           MARGIN = 1:2, 
                           FUN = mean) ## Mean Square Error
) /## Root Mean Square Error
  apply(X = cv_sys, MARGIN = 1:2, FUN = mean) ## Normalize by mean(cvs)


##################################
## Single Species STRS Optimization
##################################
target_n <- seq(from = 100, to = 40, by = -10)

index_ss_strs <- cv_ss_strs <- rb_ss_strs <- 
  array(dim = c(length(target_n), n_spp, 500),
        dimnames = list(paste0("ss_strs_n = ", target_n), spp_list, NULL))

for (ispp in 1:n_spp) { ## Loop over species -- start
  
  ## Load last run
  n_runs <- length(dir(paste0("results/", ifelse(test = igear == "otter",
                                                 yes = "otter_trawl",
                                                 no = "beam_trawl_2012_2019"),
                              "/survey_opt/SS/", spp_list[ispp], "/")))
  
  load(paste0("results/", ifelse(test = igear == "otter",
                                 yes = "otter_trawl",
                                 no = "beam_trawl_2012_2019"),
              "/survey_opt/SS/", spp_list[ispp], 
              "/Run_", n_runs, "/result_list.RData"))
  
  
  error_df <- data.frame("DOM" = "DOM1",
                         "CV1" = result_list$cvs,
                         "domainvalue"  = 1)
  
  temp_stratif <- result_list$solution$aggr_strata
  temp_stratif$DOM1 <- 1
  
  temp_bethel <- SamplingStrata::bethel(
    errors = error_df,
    stratif = temp_stratif, 
    realAllocation = T, 
    printa = T)
  temp_n <- sum(ceiling(temp_bethel))
  
  for (itarget in target_n){ ## Loop over sampling effort -- start
    while (temp_n != itarget){
      over_under <- temp_n > itarget
      CV_adj <- ifelse(over_under == TRUE, 
                       yes = 1.001,
                       no = 0.999)
      
      error_df$CV1 <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * CV_adj
      temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                            errors = error_df, 
                                            printa = TRUE)
      
      temp_n <- sum(as.numeric(temp_bethel))
      
      print(paste0("n = ", temp_n, ", CV = ", 
                   as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])) )
    }
    
    ## STRS settings
    solution <- result_list$sol_by_cell
    allocation <- as.numeric(temp_bethel)
    Nh <- tapply(X = cell_area, INDEX = solution, FUN = sum)
    Wh <- Nh / sum(Nh)
    n_strs <- sum(allocation)
    
    for (iter in 1:500) { ## Loop over iterations -- start
      
      sample_vec <- stratum_order <- c()
      for(istrata in 1:length(allocation)) { ## Loop over strata -- start
        sample_vec <- c(sample_vec,
                        sample(x = which(solution == istrata),
                               size = allocation[istrata]) )
        stratum_order <- c(stratum_order,
                           rep(x = istrata, times = allocation[istrata]))
      } ## Loop over strata -- end
      
      # Subset sub_df by which cells were chosen
      sample_df <- ms_dens[ispp, sample_vec, iter]
      
      # Calculate STRS mean density
      strata_mean <- tapply(X = sample_df, 
                            INDEX = stratum_order,
                            FUN = mean)
      tau_strs <- sum(strata_mean * Nh)
      index_ss_strs[paste0("ss_strs_n = ", itarget), ispp, iter]  <- 
        sum(strata_mean * Nh) / 1000
      
      # Calculate STRS variance of mean density
      strata_var <- tapply(X = sample_df, 
                           INDEX = stratum_order,
                           FUN = function(x)(var(x) / (length(x)) )) 
      STRS_var <- sum(Nh^2 * strata_var)
      
      # Save mean and cv of estimates across species
      cv_ss_strs[paste0("ss_strs_n = ", itarget), ispp, iter]  <- 
        sqrt(STRS_var) / tau_strs
      
      rb_ss_strs[paste0("ss_strs_n = ", itarget), , iter] <- 
        100 * (tau_strs / 1000 - true_index) / true_index
      
    } ## Loop over iterations -- end
  } ## Loop over sampling effort -- end
} ## Loop over species -- end


true_cv_ss_strs <- sweep(x = apply(X = index_ss_strs, MARGIN = 1:2, FUN = sd),
                         MARGIN = 2, STATS = true_index, FUN = "/")
rrmse_cv_ss_strs <- sqrt(apply(X = sweep(x = cv_ss_strs, 
                                         MARGIN = 1:2, 
                                         STATS = true_cv_ss_strs, 
                                         FUN = "-") ^ 2, ## Square Error
                               MARGIN = 1:2, 
                               FUN = mean) ## Mean Square Error
) /## Root Mean Square Error
  apply(X = cv_ss_strs, MARGIN = 1:2, FUN = mean) ## Normalize by mean(cvs)

##################################
## Multi Species STRS Optimization
##################################
n_runs <- length(dir(paste0("results/", ifelse(test = igear == "otter",
                                               yes = "otter_trawl",
                                               no = "beam_trawl_2012_2019"),
                            "/survey_opt/MS/")))

index_ms_strs <- cv_ms_strs <- rb_ms_strs <-
  array(dim = c(length(target_n), n_spp, 500),
        dimnames = list(paste0("ms_strs_n = ", target_n), spp_list, NULL))

load(paste0("results/", ifelse(test = igear == "otter",
                               yes = "otter_trawl",
                               no = "beam_trawl_2012_2019"),
            "/survey_opt/MS/Run_", n_runs, "/result_list.RData"))

temp_cvs <- as.numeric(result_list$cvs); names(temp_cvs) <- paste0("CV", 1:n_spp)

error_df <- cbind(data.frame("DOM" = "DOM1"),
                  as.data.frame(t(temp_cvs)),
                  data.frame("domainvalue"  = 1))

temp_stratif <- result_list$solution$aggr_strata
temp_stratif$DOM1 <- 1

temp_bethel <- SamplingStrata::bethel(
  errors = error_df,
  stratif = temp_stratif, 
  realAllocation = T, 
  printa = T)
temp_n <- sum(ceiling(temp_bethel))

for (itarget in target_n){ ## Loop over sampling effort -- start
  while (temp_n != itarget){
    over_under <- temp_n > itarget
    CV_adj <- ifelse(over_under == TRUE, 
                     yes = 1.001,
                     no = 0.999)
    
    error_df[, paste0("CV", 1:n_spp)] <- 
      as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * CV_adj
    
    temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                          errors = error_df, 
                                          printa = TRUE)
    
    temp_n <- sum(as.numeric(temp_bethel))
    
    print(paste0("n = ", temp_n, ", CV = ", 
                 as.numeric(attributes(temp_bethel)$outcv[1, "ACTUAL CV"])) )
  }
  
  ## STRS settings
  solution <- result_list$sol_by_cell
  allocation <- as.numeric(temp_bethel)
  Nh <- tapply(X = cell_area, INDEX = solution, FUN = sum)
  Wh <- Nh / sum(Nh)
  n_strs <- sum(allocation)
  
  for (iter in 1:500) { ## Loop over iterations -- start
    
    sample_vec <- stratum_order <- c()
    for(istrata in 1:length(allocation)) { ## Loop over strata -- start
      sample_vec <- c(sample_vec,
                      sample(x = which(solution == istrata),
                             size = allocation[istrata]) )
      stratum_order <- c(stratum_order,
                         rep(x = istrata, times = allocation[istrata]))
    } ## Loop over strata -- end
    
    # Subset sub_df by which cells were chosen
    sample_df <- ms_dens[, sample_vec, iter]
    
    # Calculate STRS mean density
    strata_mean <- apply(X = sample_df, 
                         MARGIN = 1, 
                         FUN = function(x) tapply(X = x, 
                                                  INDEX = stratum_order,
                                                  FUN = mean))
    
    tau_strs <- as.numeric(t(strata_mean) %*% matrix(Nh) )
    index_ms_strs[paste0("ms_strs_n = ", itarget), , iter] <- tau_strs / 1000
    
    # Calculate STRS variance of mean density
    strata_var <- apply(X = sample_df, 
                        MARGIN = 1, 
                        FUN = function(x) 
                          tapply(X = x, 
                                 INDEX = stratum_order,
                                 FUN = function(xx)
                                   (var(xx) / (length(xx)) )) )
    
    STRS_var <-  as.numeric(t(strata_var) %*% matrix(Nh^2))
    
    # Save mean and cv of estimates across species
    cv_ms_strs[paste0("ms_strs_n = ", itarget), , iter]  <- 
      sqrt(STRS_var) / tau_strs
    
    rb_ms_strs[paste0("ms_strs_n = ", itarget), , iter] <- 
      100 * (tau_strs / 1000 - true_index) / true_index
    
  } ## Loop over iterations -- end
} ## Loop over sampling effort -- end

true_cv_ms_strs <- sweep(x = apply(X = index_ms_strs, MARGIN = 1:2, FUN = sd),
                         MARGIN = 2, STATS = true_index, FUN = "/")
rrmse_cv_ms_strs <- sqrt(apply(X = sweep(x = cv_ms_strs, 
                                         MARGIN = 1:2, 
                                         STATS = true_cv_ms_strs, 
                                         FUN = "-") ^ 2, ## Square Error
                               MARGIN = 1:2, 
                               FUN = mean) ## Mean Square Error
) /## Root Mean Square Error
  apply(X = cv_ms_strs, MARGIN = 1:2, FUN = mean) ## Normalize by mean(cvs)

#############################
## Plots
#############################

{
  png(filename = "presentations/results_11_30_2021/sim_4_surveys.png", 
      width = 10, height = 10, units = "in", res = 500)
  par(mar = c(0.25, 0, 0.25, 0), mfrow = c(n_spp, 3), oma = c(4, 5, 3, 10))
  for (ispp in 1:n_spp) {
    ymax <- max(c(max(cv_srs[, ispp, ]), max(true_cv_srs[, ispp]),
                  max(cv_sys[, ispp, ]), max(true_cv_sys[, ispp]),
                  max(cv_ms_strs[, ispp, ]), max(true_cv_ms_strs[, ispp])))
    
    ## Simple Random Design
    plot(x = srs_n, y = true_cv_srs[, ispp], las = 1, type = "n",
         pch = "*", col = "red", cex = 3, ann = F, axes = F,
         ylim = c(0, ymax), 
         xlim = c(35, 105))
    box()
    axis(side = 2, las = 1)
    lines(x = srs_n, y = true_cv_srs[, ispp])
    
    if(ispp == n_spp) axis(side = 1)
    
    if (ispp == 1) mtext(side = 3, text = "Simple Random", line = 0.5, cex = 1.25)
    for (isample in 1:nrow(cv_srs)) {
      boxplot(x = cv_srs[isample, ispp, ], add = TRUE, at = srs_n[isample], 
              axes = F, width = 1, pch = 16, cex = 0.5, 
              pars = list(boxwex = 4)) 
    }
    points(x = srs_n, y = true_cv_srs[, ispp], 
           pch = "*", col = "red", cex = 3)
    lines(x = srs_n, y = true_cv_srs[, ispp])
    
    ## Systematic Design
    plot(x = sys_settings$n, y = true_cv_sys[, ispp], 
         las = 1, pch = "*", cex = 3, col = "red", axes = F, ann = F,
         ylim = c(0, ymax), xlim = c(35, 110))
    box()
    if (ispp == 1)
      mtext(side = 3, text = "Systematic",  line = 0.5, cex = 1.25)
    lines(x = sys_settings$n, y = true_cv_sys[, ispp])
    if(ispp == n_spp) axis(side = 1)
    
    for (isample in 1:nrow(sys_settings)) {
      boxplot(x = cv_sys[isample, ispp, ], 
              add = TRUE, at = sys_settings$n[isample], 
              axes = F, width = 1, pch = 16, cex = 0.5, 
              pars = list(boxwex = 4)) 
    }
    
    ## Multi Species Stratified Random Sample
    plot(x = ms_n, y = true_cv_ms_strs[, ispp], 
         las = 1, pch = "*", col = "red", axes = F, ann = F,
         ylim = c(0, ymax),  xlim = c(35, 105))
    box()
    lines(x = ms_n, y = true_cv_ms_strs[, ispp],)
    
    boxplot(t(cv_ms_strs[, ispp, ]), add = TRUE,
            at = ms_n, names = target_n, 
            pch = 16, cex = 0.25, pars = list(boxwex = 3), las = 1,
            ylim = c(0, max(cv_ms_strs[, ispp, ])),  xlim = c(35, 105), axes = F)
    points(x = ms_n, true_cv_ms_strs[, ispp], pch = "*", cex = 3, col = "red")
    if (ispp == 1) mtext(side = 3, text = "SS STRS",  line = 0.5, cex = 1.25)
    if(ispp == n_spp) axis(side = 1)
    
    
    ## Species Label
    text(x = 130, y = ymax / 2, 
         labels = gsub(x = spp_list[ispp], pattern = " ", replacement = "\n"), 
         cex = 2, xpd = NA)
    
  }
  
  mtext(side = 1, outer = TRUE, text = "Sample Size", line = 3, cex = 1.5)
  mtext(side = 2, outer = TRUE, text = "Sample CV (True CV in RED)", line = 3, cex = 1.5)
  
  dev.off()
}

#############################
## Plots: RRMSE of CV
#############################

{
  png(filename = "presentations/results_11_30_2021/sim_4_surveys_RRMSE_CV.png", 
      width = 10, height = 10, units = "in", res = 500)
  par(mar = c(0.25, 0, 0.25, 0), mfrow = c(n_spp, 3), oma = c(4, 5, 3, 10))
  for (ispp in 1:n_spp) {
    ymax <- max(c(rrmse_cv_srs[, ispp], rrmse_cv_sys[, ispp], 
                  rrmse_cv_ss_strs[, ispp], rrmse_cv_ms_strs[, ispp] )) * 1.1
    
    plot(x = srs_n, y = rrmse_cv_srs[, ispp], axes = F, 
         xlim = c(35, 105), ylim = c(0, ymax),
         pch = 16);     box()
    lines(x = srs_n, y = rrmse_cv_srs[, ispp], lwd = 1)
    if (ispp == 1) mtext(side = 3, text = "Simple Random", line = 1, cex = 1.25)
    if (ispp == n_spp) axis(side = 1)
    axis(side = 2, las = 1)
    
    plot(x = sys_settings$n, y = rrmse_cv_sys[, ispp], axes = F, ann = F,
         xlim = c(35, 105), ylim = c(0, ymax),
         pch = 16)  ; box()
    lines(x = sys_settings$n, y = rrmse_cv_sys[, ispp], lwd = 1)
    if (ispp == 1) mtext(side = 3, text = "Systematic", line = 1, cex = 1.25)
    if (ispp == n_spp) axis(side = 1)
    
    plot(x = ms_n, y = rrmse_cv_ms_strs[, ispp], axes = F, 
         xlim = c(35, 105), ylim = c(0, ymax),
         pch = 16); box()
    lines(x = ms_n, y = rrmse_cv_ms_strs[, ispp])
    if (ispp == 1) mtext(side = 3, text = "MS STRS", line = 1, cex = 1.25)
    if (ispp == n_spp) axis(side = 1)
    
    ## Species Label
    text(x = 135, y = ymax / 2,
         labels = gsub(x = spp_list[ispp], pattern = " ", replacement = "\n"),
         cex = 2, xpd = NA)
    
  }
  
  mtext(side = 1, outer = TRUE, text = "Sample Size", line = 3, cex = 1.5)
  mtext(side = 2, outer = TRUE, text = "RRMSE of CV", line = 3, cex = 1.5)
  
  dev.off()
}

goa <- sp::SpatialPointsDataFrame(coords = chukchi_grid[, c("Lon", "Lat")], proj4string = latlon_crs, data = frame_df_otter)
goa <- sp::spTransform(x = goa, CRSobj = aea_crs)
goa_ras <- raster::raster(x = goa, res = 5000)
lat_ras <- raster::rasterize(x = goa, goa_ras, field = "X1")
dist_to_shore <- raster::rasterize(x = goa, goa_ras, field = "X2")

par(mfrow = c(2, 1), mar = c(1, 1, 2, 1))
plot(lat_ras, axes = F, main = "Latitude")
plot(dist_to_shore / 1000, col = rev(magma(n = 1000)), axes = F, main = "Min. Dist. from shore (km)")

