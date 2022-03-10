###############################################################################
## Project:       Best VAST Model settings
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For each gear and species, pull the spatial and 
##                     spatiotemporal settings for the best (lowest AIC) model
###############################################################################
rm(list = ls())

for(igear in c("otter", "beam") ) { ## Loop over gear -- start
  spp_list <- dir(path = paste0("results/chukchi_", igear, "/vast_fits/"))
  
  fields_df <- data.frame()
  
  for (ispp in spp_list) {
    load(paste0("results/chukchi_", igear, "/vast_fits/", ispp, "/fit.RData"))
    field_config <- fit$data_list$FieldConfig[c("Omega", "Epsilon"), ]
    field_config <- as.vector(apply(X = field_config, 
                                    MARGIN = 1, 
                                    FUN = function(x) ifelse(x == 1, 1, 0) ))
    names(field_config) <- sapply(X = c("Spatial Field, Component ", 
                                        "Spatiotemporal Field, Component "), 
                                  FUN = function(x) paste(x, 1:2))
    
    fields_df <- rbind(fields_df, 
                       cbind(species = ispp, 
                             gear = igear,
                             t(field_config)))
  }
} ## Loop over gear -- end

