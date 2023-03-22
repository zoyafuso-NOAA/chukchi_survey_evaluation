###############################################################################
## Project:       Process NBS bottom trawl survey data for research
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), 
##                Lewis Barnett (lewis.barnett@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(dplyr)
library(RODBC)
library(readr)
library(getPass)

##################################################
####   load data from RACEBASE
##################################################
get.connected <- function(schema='AFSC'){(echo=FALSE)
  username <- getPass(msg = "Enter your ORACLE Username: ")
  password <- getPass(msg = "Enter your ORACLE Password: ")
  channel  <- RODBC::odbcDriverConnect(paste0("Driver={Oracle in OraClient12Home1};Dbq=", schema, ";Uid=", username, ";Pwd=", password, ";"))
}
channel <- get.connected()

# query directly (using function from Sean Rohan, as in coldpool package)
sql_to_rqry <- function(sql_path) {
  in_string <- readr::read_file(sql_path)
  in_string <- sub("/*.*/", "", in_string)
  out_string <- stringr::str_replace_all(in_string, pattern = "\r\n", replacement = " ") # or pattern = "\n" for some files

  return(out_string)
}

haul <- data.frame(RODBC::sqlQuery(channel, sql_to_rqry("code/data_processing/nbs_hauls_query_kotwicki_forR.sql")))
#saveRDS(haul, "data/fish_data/AK_BTS_OtterAndBeam/nbs_kotwicki_hauls.rds")
# or, once query above is conducted and saved, can load directly via
#haul <-  readRDS("data/fish_data/AK_BTS_OtterAndBeam/nbs_kotwicki_hauls.rds")

## alternative queries from Jason
# haul_jason <- RODBC::sqlQuery(channel, 
#                               sql_to_rqry("code/data_processing/NBS_EBS_query_conner.sql")) # hauls used for SAFE
# haul_jason2 <- RODBC::sqlQuery(channel, 
#                                sql_to_rqry("code/data_processing/NBS_EBS_query_conner2.sql")) # all hauls in the gear list (below)
# haul_jason2 %>% 
#   filter(START_LATITUDE > 60) %>% 
#   group_by(CRUISE) %>%
#   count()

# get cruise and catch data
# (takes a few minutes each, unless you supply a more specific catch query)
catch <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CATCH")
cruise <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CRUISE")
##################################################
####   Break down date information from haul data 
##################################################
haul$DATE <- as.Date(x = haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(x = haul$DATE)
haul$DAY <- lubridate::day(x = haul$DATE)
haul$YEAR <- lubridate::year(x = haul$DATE) 

haul %>% 
  group_by(YEAR) %>%
  count()
##################################################
####   Join with species names
##################################################
species_codes <-  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/species.csv", 
                           stringsAsFactors = FALSE)
species_codes <- dplyr::select(.data = species_codes, -YEAR_ADDED)
catch <- dplyr::left_join(x = catch, y = species_codes)
dat <- dplyr::left_join(x = haul, y = cruise, by = 'CRUISEJOIN')

##################################################
####   join with cruise data and filter NAs
##################################################
dat <- dat %>% 
  dplyr::select(-c(AUDITJOIN.x, AUDITJOIN.y, REGION.y, VESSEL.y, CRUISE.y)) %>% 
  dplyr::rename(REGION = REGION.x, 
                VESSEL = VESSEL.x, 
                CRUISE = CRUISE.x) %>% 
  dplyr::filter(!is.na(NET_WIDTH), 
                !is.na(DISTANCE_FISHED)) %>% 
  mutate(GEAR_CAT = "otter")

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
                  "Bering flounder"
                  # ,
                  # "starry flounder", 
                  # "Arctic staghorn sculpin",  "Pacific herring", 
                  # "slender eelblenny", "blue king crab", 
                  # "shorthorn (=warty) sculpin", "circumboreal toad crab", 
                  # "fuzzy hermit crab", "hairy hermit crab", "humpy shrimp", 
                  # "sculptured shrimp", "helmet crab",
                  # "Arctic argid", "kuro argid",  
                  # "green sea urchin", "notched brittlestar",
                  # "purple-orange sea star", "northern nutclam", 
                  # "common mud star", "Greenland cockle", "basketstar"
                  )

dat <- dplyr::left_join(x = dat, y = catch) %>%
  dplyr::mutate(EFFORT = DISTANCE_FISHED * NET_WIDTH * 0.001) %>%
  dplyr::mutate(CPUE_KG = WEIGHT / EFFORT, 
                CPUE_N = NUMBER_FISH / EFFORT) %>%
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
                         COMMON_NAME, SPECIES_NAME, 
                         EFFORT, CPUE_KG, CPUE_N))

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
                    c("STATION_ID", "SURVEY_NAME", "DATE", "YEAR", "MONTH", "DAY", 
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
                      "bot_depth", "bot_temp", "common_name", "catch_kg")
data_long$area_swept_km2 <- 1

##################################################
####   Save
##################################################
write.csv(x = data_long, 
          file = "data/fish_data/AK_BTS_OtterAndBeam/AK_BTS_NBS_processed_long.csv", 
          row.names = FALSE)
write.csv(x = data_wide, 
          file = "data/fish_data/AK_BTS_OtterAndBeam/AK_BTS_NBS_processed_wide.csv", 
          row.names = FALSE)

