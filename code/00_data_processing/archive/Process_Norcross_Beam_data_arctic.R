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
####  Import CPUE data
##################################################
beam_norcross <- 
  as.data.frame(read_xlsx(paste0("data/fish_data/Norcross_Beam/",
                                 "2004-09 PSBT fish data Chukchi",
                                 "_for SOAR_18Mar14.xlsx"), 
                          sheet = "CATCH 2004-2009 hauls "))

##################################################
####  Simplify field names
##################################################
beam_norcross <- subset(x = beam_norcross, 
                        select = c(HaulUniqueCDF, LatitudeStart, LongitudeStart,
                                   TempBottomC, SalBottom, DepthAvgM, Gear, 
                                   FishTaxon, Gms_per_1000sqM))
names(beam_norcross) <- c("haul_id", "lat", "lon", "bot_temp", 
                          "bot_salt", "bot_depth", "gear", "species_name", 
                          "cpue_kg_km2")

##################################################
####  Subset for species in spp_list and add common names
##################################################
beam_norcross <- subset(x = beam_norcross, subset = species_name %in% spp_list)
beam_norcross$common_name <- 
  names(spp_list)[match(beam_norcross$species_name, spp_list)]

##################################################
####  Additional manipulations:
####    add gear information, simplify names of covariates, add year
##################################################
beam_norcross$gear <- "beam"
names(beam_norcross)[4:6] <- paste0("bot_", c("temp", "salt", "depth"))

beam_norcross$year <- as.numeric(substr(x = beam_norcross$haul_id, 
                                        start = 1, stop = 4))

##################################################
####  Spread dataframe to include zeros and then lengthen the dataframe
##################################################
beam_norcross_wide <- tidyr::spread(data = subset(beam_norcross, 
                                                  select = -species_name), 
                                    value = "cpue_kg_km2", 
                                    key = "common_name",
                                    fill = 0)

beam_norcross_long <- reshape::melt(data = beam_norcross_wide, 
                                    measure.vars = names(spp_list))

names(beam_norcross_long)[9:10] <- c("common_name", "cpue_kg_km2")

##################################################
####  Save
##################################################
write.csv(x = beam_norcross_long, 
          file = "data/fish_data/Norcross_Beam/norcross_data_processed.csv", 
          row.names = FALSE)
