###############################################################################
## Project:       Synthesize Alaska Beam Trawl Data from B. Norcross
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries
##################################################
library(reshape)
library(readxl)

##################################################
####  Species to subset
##################################################
spp_list <- c("Arctic cod" = "Boreogadus saida", 
              "saffron cod" = "Eleginus gracilis",
              "walleye pollock" =  "Gadus chalcogrammus",
              "Pacific cod" = "Gadus macrocephalus",
              "Bering flounder" = "Hippoglossoides robustus",
              "yellowfin sole" = "Limanda aspera",
              "Alaska plaice" = "Pleuronectes quadrituberculatus")

##################################################
####  Import cpue data
##################################################
beam_norcross <- 
  as.data.frame(read_xlsx(paste0("data/fish_data/Norcross_Beam/",
                                 "2004-09 PSBT fish data Chukchi",
                                 "_for SOAR_18Mar14.xlsx"), 
                          sheet = "static CATCH per haul", 
                          skip = 3))

##################################################
####  Some manipulations:
####   1) some of the longitudes have the wrong sign
####   2) years as extracted out of the cruise name
##################################################
beam_norcross$lon <- ifelse(beam_norcross$LongitudeStart < 0, 
                            beam_norcross$LongitudeStart, 
                            -1 * beam_norcross$LongitudeStart) 
beam_norcross$lat <- beam_norcross$LatitudeStart

beam_norcross$year <- as.numeric(substr(x = beam_norcross$HaulUniqueCDF, 
                                        start = 1, stop = 4))

##################################################
####   "Lengthen" the dataset 
##################################################
beam_norcross_wide <- beam_norcross[, c("year", "lat", "lon", 
                                        "TempBottomC", "SalBottom", "DepthAvgM", 
                                        spp_list)]

beam_norcross_long <- reshape::melt(data = beam_norcross_wide, 
                                    measure.vars = spp_list,
                                    variable_name = "species_name")

##################################################
####  Additional manipulations: add common names and a field to denote that 
####  this is beam data
##################################################
beam_norcross_long$common_name <- 
  names(spp_list)[match(beam_norcross_long$species_name, spp_list)]
beam_norcross_long$gear <- "beam"
names(beam_norcross_long)[c(4:6, 8)] <- 
  c(paste0("bot_", c("temp", "salt", "depth")), "catch_kg")

##################################################
####  Save
##################################################
write.csv(x = beam_norcross_long, 
          file = "data/fish_data/Norcross_Beam/norcross_data_processed.csv", 
          row.names = FALSE)
