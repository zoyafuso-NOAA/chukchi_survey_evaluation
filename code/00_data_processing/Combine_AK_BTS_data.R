###############################################################################
## Topic:         Combine bottom trawl survey data for NBS-Chukchi research
## Author:        Lewis Barnett (lewis.barnett@noaa.gov)
###############################################################################

library(tidyverse)

list_of_files <- c(list.files(path = "data/fish_data/AK_BTS_OtterAndBeam/data_long_by_taxa/",
                            pattern = "\\.csv$",
                            full.names = TRUE),
                   list.files(path = "data/fish_data/AK_BTS_OtterAndBeam/data_long_by_taxa_nbs/",
                              pattern = "\\.csv$",
                              full.names = TRUE))

df <- list_of_files %>% 
  map_dfr(read.csv, header=TRUE, fill=TRUE)