###############################################################################
## Project:       Simulate Surveys Function
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create function to calcualte a STRS and output mean, 
##                CV, and relative bias
###############################################################################
rm(list = ls())

library(sp)
library(rgeos)

##################################
## Simulate systematic survey based on 2012 otter Trawl Design
##################################
load("data/survey_opt_data/optimization_data_otter.RData")

chukchi_grid <- read.csv(paste0("data/spatial_data/",
                                "BS_Chukchi_extrapolation_grids/",
                                "ChukchiThorsonGrid.csv"))
chukchi_grid$Area_km2 <- chukchi_grid$Shape_Area / 1000 / 1000

otter <- subset(x = read.csv(paste0("data/fish_data/AK_BTS_OtterAndBeam/",
                                    "AK_BTS_Arctic_processed_wide.csv")),
                subset = GEAR_CAT == "otter" & YEAR == 2012)
otter_coords <- 
  sp::spTransform(x = sp::SpatialPoints(coords = otter[, c("MEAN_LONGITUDE", 
                                                           "MEAN_LATITUDE")], 
                                        proj4string = latlon_crs), 
                  CRSobj = aea_crs) 

which_grid_pts <- apply(X = rgeos::gDistance(spgeom1 = otter_coords, 
                                             spgeom2 = grid_pts, 
                                             byid = TRUE), 
                        MARGIN = 2, 
                        FUN = which.min)

spp_list
ispp = 4
n_runs <- length(dir(paste0("results/otter_trawl/survey_opt/SS/", 
                            spp_list[ispp], "/")))
cv_n <- data.frame(run = 1:n_runs)

for (irun in 1:n_runs){
  load(paste0("results/otter_trawl/survey_opt/SS/", spp_list[ispp], 
              "/Run_", irun, "/result_list.RData"))
  cv_n$cv[irun] <- result_list$cvs
  cv_n$n[irun] <- result_list$n
}

idx <- which.min(abs(cv_n$n - 71))
cv_n[idx, ]

load(paste0("results/otter_trawl/survey_opt/SS/", spp_list[ispp], 
            "/Run_", idx, "/result_list.RData"))

load(paste0("results/otter_trawl/", spp_list[ispp], "/simulated_data.RData"))

str(sim_data)

n = sum(chukchi_grid$Area_km2[which_grid_pts])
N = sum(chukchi_grid$Area_km2)
srs_cv <- strs_cv <- sys_cv <- srs_index <- strs_index <- sys_index <- c()

Nh <- result_list$sum_stats$Population
Wh <- Nh / sum(Nh)
cell_area <- chukchi_grid$Shape_Area / 1000 / 1000
solution <- result_list$sol_by_cell
nh <- result_list$sum_stats$Allocation
wh <- nh / Nh

for (iter in 1:500) {
  
  # Simple Random Design
  set.seed(iter * 234)
  srs_sample <- sample(x = sim_data[, iter], size = length(otter_coords))
  srs_tau <- N * mean(srs_sample)
  srs_sd <- sqrt(N * (N - n) * var(srs_sample) / n)
  
  srs_index[iter] <- srs_tau / 1000
  srs_cv[iter] <- srs_sd / srs_tau
  
  # Systematic Design
  Y <- sim_data[which_grid_pts, iter]
  sys_sd <- sqrt(N * (N - n) * var(Y) / n)
  sys_cv[iter] <- sys_sd / ( mean(Y) * N)
  sys_index[iter] <-  mean(Y) * N / 1000
  
  # Take a random sample based on the allocation and stratum
  sample_vec <- c()
  for(istrata in 1:length(nh)) {
    sample_vec <- c(sample_vec,
                    sample(x = which(solution == istrata),
                           size = nh[istrata]) )
  }
  
  sampled_strata <- rep(x = 1:length(nh), 
                        times = nh)
  
  # Subset sub_df by which cells were chosen
  sample_df <- sim_data[sample_vec, iter]
  
  # Calculate STRS mean density
  strata_mean <- tapply(X = sample_df, 
                        INDEX = sampled_strata,
                        FUN = mean)
  STRS_mean <- sum(strata_mean * Wh)
  
  # Calculate STRS variance of mean density
  strata_var <- tapply(X = sample_df, 
                       INDEX = sampled_strata,
                       FUN = var)
  STRS_var <- sum(strata_var * Wh^2 * (1 - wh) / nh)
  
  # Save mean and cv of estimates across species
  strs_cv[iter] <- sqrt(STRS_var) / STRS_mean
  
  strs_index[iter] <- sum(strata_mean * tapply(cell_area, solution, sum)) / 1000
}

boxplot(cbind((srs_index - true_index[ispp]) / true_index[ispp], 
              (sys_index - true_index[ispp]) / true_index[ispp],
              (strs_index - true_index[ispp]) / true_index[ispp]))
abline(h = 0)

true_cv <- apply(X = cbind(srs_index, sys_index, strs_index), 
                 MARGIN = 2, 
                 FUN = function(x) 
                   sd(x - true_index[ispp]) / true_index[ispp] )

rrmse_cv <- sweep(x = cbind(srs_cv, sys_cv, strs_cv), 
                  MARGIN = 2, 
                  FUN = "-", 
                  STATS = true_cv)
rrmse_cv <- apply(X = rrmse_cv, 
                  MARGIN = 2, 
                  FUN = function(x) sqrt(mean(x^2)) )
rrmse_cv <- rrmse_cv / colMeans(cbind(srs_cv, sys_cv, strs_cv))


boxplot(cbind(srs_cv, sys_cv, strs_cv))
abline(h = true_cv, lwd= 2, col = c("red", "black", "blue"))
rrmse_cv
