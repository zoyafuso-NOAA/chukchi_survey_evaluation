##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Best VAST Model settings
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For each gear and species, pull the spatial and 
##                     spatiotemporal settings for the best (lowest AIC) model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

fields_df <- 
  rbind(
    data.frame("Taxon" = dir(path = "results/chukchi_otter/vast_fits/"),
               "Gear" = "otter"),
    data.frame("Taxon" = dir(path = "results/chukchi_beam/vast_fits/"),
               "Gear" = "beam"))

fields_df <- fields_df[order(fields_df$Taxon), ]
fields_df[, sapply(X = c("L_omega", 
                         "L_epsilon"), 
                   FUN = function(x) paste0(x, 1:2, "_z"))] <- NA


for(irow in 1:nrow(fields_df) ) { ## Loop over gear -- start
  
  vast_fit <- paste0("results/chukchi_", fields_df$Gear[irow], 
                     "/vast_fits/", fields_df$Taxon[irow], "/fit.RData")
  if (file.exists(vast_fit)) {
    load(vast_fit)
    
    for (ifield in c("omega", "epsilon")) {
      for (icomponent in 1:2) {
        fields_df[irow, paste0("L_", ifield, icomponent, "_z")] <- 
          ifelse(test = paste0("L_", ifield, icomponent, "_z") %in% 
                   names(fit$parameter_estimates$par),
                 yes = "x", no = "")
      }
    }
  }
} ## Loop over gear -- end

