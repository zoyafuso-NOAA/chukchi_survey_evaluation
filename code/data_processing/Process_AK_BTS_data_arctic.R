###############################################################################
## Project:       Synthesize Alaska Bottom Trawl Arctic otter trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), modfied from 
##                Lewis Barnett (lewis.barnett@noaa.gov)
##
## Notes:         2012 is only year where both beam and 83-112 otter trawls 
##                were conducted in the Chukchi. 83-112 trawl data includes 
##                many values of GEAR: c(30:44,47,178,180,181). Code 30 from 
##                1975-1982/1985 then small changes through to code 44 used for
##                current EBS survey; code 47 not present in database. Code 178
##                and 180 from 1998 only; code 181 not present in database. 
##                Beam trawl is GEAR_CODE 180 (when fished as extension of 
##                otter trawl, unclear whether code is for otter trawl only or
##                not), all from year 1998 in BS REGION (also 179 in codebook 
##                but not present in data). Plum staff beam trawl primarily 218 
##                and 219 (Norcross; but only have year 2012, when also used 
##                gear code 44), 217 also similar (but only have handful in 
##                1989 for BS region)
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(dplyr)

##################################################
####   load flat files extracted from RACEBASE
##################################################
cruise <- read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/cruise.csv",
                   stringsAsFactors = FALSE)
haul <- read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/haul.csv",
                 stringsAsFactors = FALSE)
catch <- read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/catch_subsetted_BS.csv", 
                  stringsAsFactors = FALSE)

##################################################
####   Some manipulations of the haul dataframe:
####
####   Break down date information from haul data (note that year conversion 
####   for years prior to 1969 are incorrectly assigned to 2000s)
####   
####   Filter to include HAUL_TYPE 3 (standard planned haul), 0 (opportunistic,
####   haul), and 23 (gear performance hauls)
####
####   Include Good Performance hauls (PERFORMANCE >= 0)
##################################################
haul$DATE <- as.Date(x = haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(x = haul$DATE)
haul$DAY <- lubridate::day(x = haul$DATE)
haul$YEAR <- lubridate::year(x = haul$DATE) 

haul <- subset(x = haul, 
               subset = PERFORMANCE >= 0 
               & HAUL_TYPE %in% c(0, 3, 23))

##################################################
####   Join with species names
##################################################
species_codes <- 
  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/species.csv", 
           stringsAsFactors = FALSE)
species_codes <- dplyr::select(.data = species_codes, -YEAR_ADDED)
catch <- dplyr::left_join(x = catch, y = species_codes)
dat <- dplyr::left_join(x = haul, y = cruise, by = 'CRUISEJOIN')

##################################################
####   Filter by gear, region, latitude
##################################################
dat <- dat %>% 
  dplyr::select(-c(AUDITJOIN.x, AUDITJOIN.y, REGION.y, VESSEL.y, CRUISE.y)) %>% 
  dplyr::rename(REGION = REGION.x, 
                VESSEL = VESSEL.x, 
                CRUISE = CRUISE.x) %>% 
  dplyr::filter(REGION == "BS", 
                GEAR %in% c(30:44, 47, 178, 180, 181, 217, 218, 219), 
                START_LATITUDE > 65.5 | END_LATITUDE > 65.5) %>% 
  mutate(GEAR_CAT = dplyr::case_when(GEAR %in% c(30:44, 47, 
                                                 178, 180, 181) ~ "otter", 
                                     GEAR %in% c(179, 217, 
                                                 218, 219)~ "beam"))

dat$GEAR_CAT <- as.factor(dat$GEAR_CAT)

##################################################
####   Join haul and catch data for selected species in species_list, 
####   summarize CPUE in numbers or kg per square km currently 
##################################################
species_list <- c("Arctic cod", "saffron cod", "Pacific cod", "walleye pollock",
                  "snow crab", "yellowfin sole","Alaska plaice", 
                  "Bering flounder")

dat <- dplyr::left_join(x = dat, y = catch) %>%
  dplyr::mutate(CPUE_KG = WEIGHT / (DISTANCE_FISHED * NET_WIDTH * 0.001), 
                CPUE_N = NUMBER_FISH / (DISTANCE_FISHED * NET_WIDTH * 0.001)) %>%
  dplyr::filter(COMMON_NAME %in% species_list,
                !is.na(CPUE_KG),
                !is.na(CPUE_N))

##################################################
#### Calculate station location as the mean of the start and ending points
##################################################
dat$MEAN_LATITUDE <- rowMeans(dat[, c("START_LATITUDE", "END_LATITUDE")])
dat$MEAN_LONGITUDE <- rowMeans(dat[, c("START_LONGITUDE", "END_LONGITUDE")])

##################################################
#### Subset gear and year-specific data
##################################################
dat <- subset(dat,
              select = c(SURVEY_NAME, YEAR, GEAR_CAT,
                         STATIONID, MEAN_LONGITUDE, MEAN_LATITUDE,
                         DATE, MONTH, DAY, YEAR,
                         GEAR_DEPTH, GEAR_TEMPERATURE,
                         COMMON_NAME, SPECIES_NAME, CPUE_KG, CPUE_N))

##################################################
#### Create a wide df with zeros filled for stations where species not observed
##################################################
dat$STATION_ID <- with(dat, paste(MEAN_LONGITUDE, MEAN_LATITUDE))

cpue_wide <- tidyr::spread(
  data = dat[, c("YEAR", "MONTH", "DAY", "DATE",
                 "STATION_ID", "CPUE_KG", "COMMON_NAME", "GEAR_CAT")], 
  key = COMMON_NAME, 
  value = "CPUE_KG", 
  fill = 0)

cpue_wide <- cpue_wide[order(cpue_wide$GEAR_CAT,
                             cpue_wide$STATION_ID), ]

##################################################
## attach station specific data to data_wide
##################################################
station_data <- dat[!duplicated(dat$STATION_ID), 
                    c("STATION_ID", "SURVEY_NAME", 
                      "DATE", "YEAR", "MONTH", "DAY", 
                      "GEAR_CAT", "GEAR_DEPTH", "GEAR_TEMPERATURE",  
                      "MEAN_LONGITUDE", "MEAN_LATITUDE") ]

station_data <- station_data[order(station_data$GEAR_CAT,
                                   station_data$STATION_ID), ]

data_wide <- cbind(station_data, 
                   subset(cpue_wide, select = -STATION_ID))

##################################################
#### "lengthen" wide dataset to long-form to match VAST data input 
##################################################
data_long <- reshape::melt(
  data = data_wide[, c("YEAR", "MONTH", "DAY", "DATE", "GEAR_CAT", 
                       "MEAN_LONGITUDE", "MEAN_LATITUDE", 
                       "GEAR_DEPTH", "GEAR_TEMPERATURE",
                       unique(dat$COMMON_NAME))],
  measure.vars = unique(dat$COMMON_NAME),
  variable_name = "COMMON_NAME")

names(data_long) <- c("year", "month", "day", "date", "gear", "lon", "lat", 
                      "bot_depth", "bot_temp", "common_name", "cpue_kg_km2")
data_long$area_swept_km2 <- 1

##################################################
####   Save
##################################################
write.csv(x = data_long, 
          file = "data/fish_data/AK_BTS_OtterAndBeam/AK_BTS_Arctic_processed_long.csv", 
          row.names = FALSE)
write.csv(x = data_wide, 
          file = "data/fish_data/AK_BTS_OtterAndBeam/AK_BTS_Arctic_processed_wide.csv", 
          row.names = FALSE)

