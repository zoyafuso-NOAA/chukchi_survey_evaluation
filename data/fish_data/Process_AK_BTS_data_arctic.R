###############################################################################
## Project:       Get AK survey data for arctic areas
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
####   Set up directories
##################################################
wd <- "C:/Users/zack.oyafuso/Desktop/data-raw/"
output_wd <- "C:/Users/zack.oyafuso/Desktop/Arctic/data/"

##################################################
####   Import Libraries
##################################################
library(dplyr)

##################################################
####   load flat files extracted from RACEBASE
##################################################
cruise <- read.csv(file = paste0(wd, "cruise.csv"), stringsAsFactors = FALSE)
haul <- read.csv(file = paste0(wd, "haul.csv"), stringsAsFactors = FALSE)
catch <- read.csv(file = paste0(wd, "catch.csv"), stringsAsFactors = FALSE)

##################################################
####   Break down date information from haul data (note that year conversion 
####   for years prior to 1969 are incorrectly assigned to 2000s)
##################################################
haul$DATE <- as.Date(x = haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(x = haul$DATE)
haul$DAY <- lubridate::day(x = haul$DATE)
haul$YEAR <- lubridate::year(x = haul$DATE) 

##################################################
####   Join with species names
##################################################
species_codes <-  read.csv(file = paste0(wd, "species.csv"), 
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

# exploring catches by species
# catch_rank <- dplyr::left_join(dat, catch) %>%
#   dplyr::group_by(COMMON_NAME) %>%
#   dplyr::summarize(sum_n = sum(NUMBER_FISH, na.rm = TRUE)) %>%
#   dplyr::arrange(desc(sum_n))
# print(catch_rank, n=300)

##################################################
####   Join haul and catch data for selected species, summarize CPUE in numbers
####   or kg per square km currently including shrimps, crabs, and other mobile 
####   (non-molluscan) epifauna in top 40 total catch numbers, along with blue 
####   king crab
##################################################
species_list <- c("Arctic cod", "saffron cod", "Pacific cod", "walleye pollock",
                  "snow crab", "yellowfin sole","Alaska plaice", 
                  "Bering flounder", "starry flounder", 
                  "Arctic staghorn sculpin",  "Pacific herring", 
                  "slender eelblenny", "blue king crab", 
                  "shorthorn (=warty) sculpin", "circumboreal toad crab", 
                  "fuzzy hermit crab", "hairy hermit crab", "humpy shrimp", 
                  "sculptured shrimp", "helmet crab",
                  "Arctic argid", "kuro argid",  
                  "green sea urchin", "notched brittlestar",
                  "purple-orange sea star", "northern nutclam", 
                  "common mud star", "Greenland cockle", "basketstar")

# c("rainbow smelt", "circumpolar eualid", "shortscale eualid", 
# "smooth nutclam", "brownscaled sea cucumber")

dat <- dplyr::left_join(x = dat, y = catch) %>%
  dplyr::mutate(CPUE_KG = WEIGHT / (NET_HEIGHT * NET_WIDTH * 0.001), 
                CPUE_N = NUMBER_FISH / (NET_HEIGHT * NET_WIDTH * 0.001)) %>%
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
                         GEAR_DEPTH, GEAR_TEMPERATURE,
                         COMMON_NAME, SPECIES_NAME, CPUE_KG, CPUE_N))

##################################################
#### Create a wide df with zeros filled for stations where species not observed
##################################################
dat$STATION_ID <- with(dat, paste(MEAN_LONGITUDE, MEAN_LATITUDE))

cpue_wide <- tidyr::spread(
  data = dat[, c("YEAR", "STATION_ID", "CPUE_KG", "COMMON_NAME", "GEAR_CAT")], 
  key = COMMON_NAME, 
  value = "CPUE_KG", 
  fill = 0)

cpue_wide <- cpue_wide[order(cpue_wide$GEAR_CAT,
                             cpue_wide$STATION_ID), ]

##################################################
## attach station specific data to data_wide
##################################################
station_data <- dat[!duplicated(dat$STATION_ID), 
                    c("SURVEY_NAME", "YEAR", "GEAR_CAT", "STATION_ID", 
                      "MEAN_LONGITUDE", "MEAN_LATITUDE") ]
station_data <- station_data[order(station_data$GEAR_CAT,
                                   station_data$STATION_ID), ]

data_wide <- cbind(station_data, 
                   subset(cpue_wide, select = -STATION_ID))

##################################################
#### "lengthen" wide dataset to long-form to match VAST data input 
##################################################
data_long <- reshape::melt(
  data = data_wide[, c("YEAR", "GEAR_CAT", "MEAN_LONGITUDE", 
                       "MEAN_LATITUDE", unique(dat$COMMON_NAME))],
  measure.vars = unique(dat$COMMON_NAME),
  variable_name = "COMMON_NAME")

names(data_long) <- c("year", "gear", "lon", "lat", "common_name", "catch_kg")
data_long$area_swept_km2 <- 1

##################################################
####   Save
##################################################
write.csv(x = data_long, 
          file = paste0(output_wd, "AK_BTS_ARCTIC.csv"), 
          row.names = FALSE)
