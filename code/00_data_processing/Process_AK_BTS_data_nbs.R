###############################################################################
## Project:       Process NBS 83-112 bottom trawl survey data for research
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

sql_to_rqry <- function(sql_path) {
  in_string <- readr::read_file(sql_path)
  in_string <- sub("/*.*/", "", in_string)
  out_string <- stringr::str_replace_all(in_string, pattern = "\r\n", replacement = " ") # or pattern = "\n" for some files
  return(out_string)
}

haul <- data.frame(RODBC::sqlQuery(channel, sql_to_rqry("code/00_data_processing/nbs_hauls_query_kotwicki_forR.sql")))
#saveRDS(haul, "data/fish_data/AK_BTS_OtterAndBeam/nbs_kotwicki_hauls.rds")
# or, once query above is conducted and saved, can load directly via
#haul <-  readRDS("data/fish_data/AK_BTS_OtterAndBeam/nbs_kotwicki_hauls.rds")

# get cruise and catch data
# (takes a few minutes each, unless you supply a more specific catch query)
catch <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CATCH")
cruise <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CRUISE")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create result directory ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
output_dir <- "data/fish_data/AK_BTS_OtterAndBeam/data_long_by_taxa_nbs/"
if (!dir.exists(paths = output_dir)) dir.create(path = output_dir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Manipulate haul data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
haul$DATE <- as.Date(x = haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(x = haul$DATE)
haul$DAY <- lubridate::day(x = haul$DATE)
haul$YEAR <- lubridate::year(x = haul$DATE) 

haul %>% 
  group_by(YEAR) %>%
  count()

# filter to YEAR > 2000 for A-IERP? 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Join haul data with species names ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
species_codes <- 
  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/species.csv", 
           stringsAsFactors = FALSE)
taxa <- read.csv(
  file = "data/fish_data/AK_BTS_OtterAndBeam/taxonomy.csv", 
  stringsAsFactors = FALSE)
species_codes <- dplyr::select(.data = species_codes, -YEAR_ADDED)
catch <- dplyr::left_join(x = catch, y = species_codes)
dat <- dplyr::left_join(x = haul, y = cruise, by = 'CRUISEJOIN')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Combine haul and cruise data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
dat <- dat %>% 
  dplyr::select(-c(AUDITJOIN.x, AUDITJOIN.y, REGION.y, VESSEL.y, CRUISE.y)) %>% 
  dplyr::rename(REGION = REGION.x, 
                VESSEL = VESSEL.x, 
                CRUISE = CRUISE.x)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate CPUE and effort ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
dat <- dplyr::left_join(x = dat, y = catch) %>%
  dplyr::mutate(CPUE_KG = WEIGHT / (DISTANCE_FISHED * NET_WIDTH * 0.001),
                CPUE_N = NUMBER_FISH / (DISTANCE_FISHED * NET_WIDTH * 0.001))

dat$AREA_SWEPT <- dat$DISTANCE_FISHED * dat$NET_WIDTH * 0.001
dat <- dat[!is.na(x = dat$AREA_SWEPT), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Station centroids ----
##   Calculate station location as the mean of the start and ending points
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
dat$MEAN_LATITUDE <- rowMeans(x = dat[, c("START_LATITUDE", "END_LATITUDE")])
dat$MEAN_LONGITUDE <- rowMeans(x = dat[, c("START_LONGITUDE", "END_LONGITUDE")])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Subset year-specific data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
dat <- subset(x = dat,
              select = c(SURVEY_NAME, YEAR, HAULJOIN,
                         STATIONID, MEAN_LONGITUDE, MEAN_LATITUDE,
                         DATE, MONTH, DAY, YEAR,
                         GEAR_DEPTH, GEAR_TEMPERATURE,
                         SPECIES_CODE, 
                         CPUE_KG, CPUE_N, AREA_SWEPT))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Add zeros to df ----
##   Create a wide df with zeros filled for stations where species not observed
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_wide <- tidyr::spread(
  data = dat[, c("YEAR", "MONTH", "DAY", "DATE",
                 "HAULJOIN", "CPUE_KG", "AREA_SWEPT",
                 "SPECIES_CODE")], 
  key = SPECIES_CODE, 
  value = "CPUE_KG", 
  fill = 0)

cpue_wide <- cpue_wide[order(cpue_wide$HAULJOIN), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Attach station specific data to data_wide ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
station_data <- dat[!duplicated(dat$HAULJOIN) ,
                    c("HAULJOIN", "SURVEY_NAME", 
                      "DATE", "YEAR", "MONTH", "DAY", 
                      "GEAR_DEPTH", "GEAR_TEMPERATURE",  
                      "MEAN_LONGITUDE", "MEAN_LATITUDE") ]

station_data <- station_data[order(station_data$HAULJOIN), ]

data_wide <- cbind(station_data, cpue_wide)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species-specific data outut ----
##   Loop over individual species, save data inputs to file and remove from
##      from the cpue_wide. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_info <- data.frame()

solo_spp <- subset(x = species_codes, 
                   subset = COMMON_NAME %in% 
                     c("Arctic cod",  "Pacific cod", "walleye pollock", 
                       "snow crab", "yellowfin sole", "Alaska plaice", 
                       "Bering flounder", "saffron cod",
                       "Pacific herring",
                       "purple-orange sea star") )

for (irow in 1:nrow(solo_spp)) {
  
  ispp_code <- solo_spp$SPECIES_CODE[irow]
  ispp_sci_name <- solo_spp$SPECIES_NAME[irow]
  ispp <- solo_spp$COMMON_NAME[irow]
  
  species_info <- rbind(species_info, 
                        data.frame(common_name = ispp,
                                   scientific_name = ispp_sci_name,
                                   file_name = ispp,
                                   species_code = ispp_code))
  
  if (ispp_code %in% names(data_wide)) {
    data_long <- data_wide[, c("YEAR", "MONTH", "DAY", "DATE", 
                               "MEAN_LONGITUDE", "MEAN_LATITUDE", "HAULJOIN", 
                               "GEAR_DEPTH", "GEAR_TEMPERATURE", 
                               "AREA_SWEPT", ispp_code)]
    names(data_long) <- c("year", "month", "day", "date", "lon", "lat", 
                          "hauljoin", "bot_depth", "bot_temp", 
                          "actual_area_swept", "cpue_kg_km2")
    
    ## Remove species column from data_wide
    data_wide <- data_wide[!names(data_wide) %in% paste(ispp_code)]
    
    ## Save data
    write.csv(x = data_long, file = paste0(output_dir, ispp, ".csv"), 
              row.names = F) 
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Aggregate Species ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
aggregate_species <- read.csv(file = "data/fish_data/aggregates_species.csv")
for (irow in 1:nrow(aggregate_species)) {
  
  species_code_query <- paste0("subset(x = taxa, subset = ", 
                               aggregate_species$species_query[irow], ")")
  species_query <- eval(parse(text = species_code_query))
  species_query <- 
    species_query[species_query$SPECIES_CODE %in% names(data_wide), ]
  
  if (nrow(species_query) > 0) {
    ispp <- aggregate_species$taxon_name[irow]
    ispp_code <- species_query$SPECIES_CODE
    ispp_common <- species_query$REPORT_NAME_SCIENTIFIC 
    ispp_sci_name <- species_query$REPORT_NAME_SCIENTIFIC
    
    species_info <- rbind(species_info, 
                          data.frame(common_name = ispp,
                                     scientific_name = ispp_sci_name,
                                     file_name = ispp,
                                     species_code = ispp_code))
    
    
    sub_df <- as.matrix(data_wide[, paste(ispp_code)])
    
    data_long <- cbind(data_wide[, c("YEAR", "MONTH", "DAY", "DATE", 
                                     "MEAN_LONGITUDE", 
                                     "MEAN_LATITUDE", "HAULJOIN", 
                                     "GEAR_DEPTH", "GEAR_TEMPERATURE", 
                                     "AREA_SWEPT")],
                       cpue_kg_km2 = rowSums(sub_df))
    names(data_long) <- c("year", "month", "day", "date", "lon", "lat",
                          "hauljoin", "bot_depth", "bot_temp", 
                          "actual_area_swept", "cpue_kg_km2")
    
    ## Remove species column from data_wide
    data_wide <- data_wide[!names(data_wide) %in% paste(ispp_code)]
    
    ## Save data
    write.csv(x = data_long, file = paste0(output_dir, ispp, ".csv"),
              row.names = F) 
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save species info data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(x = species_info, 
          file = paste0(dirname(output_dir), "/species_info.csv"),
          row.names = F)