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
fish_len <- read.csv(file = "C:/Users/zack.oyafuso/Desktop/data-raw/length.csv")

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
fish_len <- dplyr::left_join(x = fish_len, y = species_codes)
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

otter_lens <- dplyr::left_join(x = dat, y = fish_len) %>%
  dplyr::filter(COMMON_NAME %in% species_list)

##################################################
#### Calculate station location as the mean of the start and ending points
##################################################
otter_lens$MEAN_LATITUDE <- rowMeans(otter_lens[, c("START_LATITUDE", 
                                                    "END_LATITUDE")])
otter_lens$MEAN_LONGITUDE <- rowMeans(otter_lens[, c("START_LONGITUDE",
                                                     "END_LONGITUDE")])

##################################################
#### Save
##################################################
write.csv(x = otter_lens, 
          file = "data/fish_data/AK_BTS_OtterAndBeam/AK_BTS_Arctic_processed_lengths_long.csv",
          row.names = FALSE)

