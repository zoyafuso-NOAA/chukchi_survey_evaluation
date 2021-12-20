###############################################################################
## Project:       Synthesize Alaska Bottom Trawl Arctic otter trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), modfied from 
##                Lewis Barnett (lewis.barnett@noaa.gov)
##
## Notes:         2008 Beaufort Otter Trawl Data
##                Use 32-bit of R
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(RODBC)
library(tidyr)

##################################################
####   Import Data
##################################################
data_channel <- RODBC::odbcConnectAccess("data/fish_data/2008_Beaufort/BSS_Database2008_Correct.mdb")

station_locs <- sqlFetch(channel = data_channel, 
                         sqtable = "tblAllStationLocations")
catch <- sqlFetch(channel = data_channel, 
                  sqtable = "tblBotCatch")
haul <- sqlFetch(channel = data_channel, 
                 sqtable = "tblBotHaul")

close(data_channel)

##################################################
####   Species List
##################################################
spp_list <- c("Arctic cod" = "Boreogadus saida", 
              # "saffron cod" = "Eleginus gracilis",
              # "walleye pollock" =  "Gadus chalcogrammus",
              "walleye pollock" = "Theragra chalcogramma",
              # "Pacific cod" = "Gadus macrocephalus",
              "Bering flounder" = "Hippoglossoides robustus"#,
              # "yellowfin sole" = "Limanda aspera",
              # "Alaska plaice" = "Pleuronectes quadrituberculatus"
)

##################################################
####   Data Manipulations:
####    1) subset haul df to only bottom trawls with satisfactory 
####       performance (performance >= 0)
####    2) based on 1, subset station_locs and catch dfs to only thoses 
####          with haul numbers in haul and species in the spp_list
####    3) attach area_swept to haul locations
##################################################
haul <- subset(x = haul, 
               subset = Performance >= 0)
station_locs <- subset(x = station_locs, 
                       subset = EventDescription == "Bottom" &
                         Haul %in% haul$Haul)
catch <- subset(catch, 
                subset = Haul %in% haul$Haul &
                  SpeciesLatinName %in% spp_list,
                select = c(Haul, SpeciesLatinName, TotalWeight))

haul[, c("Lon", "Lat")] <- station_locs[match(haul$Haul, station_locs$Haul), 
                                        c("LongDD", "LatDD")]

##################################################
####   Widen catch df to include zeros and then relengthen df
##################################################
data_wide <- cbind(Year = substr(x = haul$Cruise, start = 1, stop = 4),
                   subset(x = haul, 
                          select = c(NetType, AreaSweptkm2, 
                                     Lon, Lat, AvgBottomDepth)),
                   
                   tidyr::spread(data = catch, 
                                 value = TotalWeight, 
                                 key = SpeciesLatinName, 
                                 fill = 0))

data_long <- reshape::melt(data = data_wide, 
                           id.vars = 1:7, 
                           variable_name = "species_name")
names(data_long) <- c("year", "net_type", "area_swept_km2", "lon", "lat", 
                      "avg_bottom_depth", "haul", "species_name", "biomass_kg")

##################################################
####   Add common names and cpue calculation
##################################################
data_long$common_name <- 
  names(spp_list)[match(data_long$species_name, spp_list)]
data_long$cpue_kg_km2 <- data_long$biomass_kg / data_long$area_swept_km2
data_long$gear <- "otter"

##################################################
####   Save
##################################################
write.csv(x = data_long,
          file = paste0("data/fish_data/2008_Beaufort/",
                        "AK_BTS_Beaufort_processed_long.csv"),
          row.names = FALSE)
