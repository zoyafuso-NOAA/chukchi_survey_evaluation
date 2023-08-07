##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Synthesize Alaska Bottom Trawl Arctic otter trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Lewis Barnett (lewis.barnett@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(dplyr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load flat files extracted from RACEBASE ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
cruise <- 
  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/cruise.csv",
           stringsAsFactors = FALSE)
haul <- 
  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/haul.csv",
           stringsAsFactors = FALSE)
catch <- 
  read.csv(file = "data/fish_data/AK_BTS_OtterAndBeam/catch_subsetted_BS.csv", 
           stringsAsFactors = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create result directory ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
output_dir <- "data/fish_data/AK_BTS_OtterAndBeam/data_long_by_taxa/"
if (!dir.exists(paths = output_dir)) dir.create(path = output_dir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Manipulate haul data ----
##
##   Break down date information from haul data (note that year conversion 
##   for years prior to 1969 are incorrectly assigned to 2000s)
##   
##   Filter to include HAUL_TYPE 3 (standard planned haul), 0 (opportunistic,
##   haul), and 23 (gear performance hauls)
##
##   Include Good Performance hauls (PERFORMANCE >= 0)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
haul$DATE <- as.Date(x = haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(x = haul$DATE)
haul$DAY <- lubridate::day(x = haul$DATE)
haul$YEAR <- lubridate::year(x = haul$DATE) 

haul <- subset(x = haul, 
               subset = PERFORMANCE >= 0 & HAUL_TYPE %in% c(0, 3, 23))

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
##   Filter by gear, region, latitude ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
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
##   Subset gear and year-specific data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
dat <- subset(x = dat,
              select = c(SURVEY_NAME, YEAR, GEAR_CAT, HAULJOIN,
                         STATIONID, MEAN_LONGITUDE, MEAN_LATITUDE,
                         DATE, MONTH, DAY, YEAR,
                         GEAR_DEPTH, GEAR_TEMPERATURE,
                         SPECIES_CODE, WEIGHT,
                         CPUE_KG, CPUE_N, AREA_SWEPT))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Add zeros to df ----
##   Create a wide df with zeros filled for stations where species not observed
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_wide <- tidyr::spread(
  data = dat[, c("YEAR", "MONTH", "DAY", "DATE",
                 "HAULJOIN", "WEIGHT", "CPUE_KG", "CPUE_N","AREA_SWEPT",
                 "SPECIES_CODE", "GEAR_CAT")], 
  key = SPECIES_CODE, 
  value = "CPUE_KG", 
  fill = 0)

cpue_wide <- cpue_wide[order(cpue_wide$GEAR_CAT,
                             cpue_wide$HAULJOIN), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Attach station specific data to data_wide ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
station_data <- dat[!duplicated(dat$HAULJOIN) ,
                    c("HAULJOIN", "SURVEY_NAME", 
                      "DATE", "YEAR", "MONTH", "DAY", 
                      "GEAR_CAT", "GEAR_DEPTH", "GEAR_TEMPERATURE",  
                      "MEAN_LONGITUDE", "MEAN_LATITUDE") ]

station_data <- station_data[order(station_data$GEAR_CAT,
                                   station_data$HAULJOIN), ]

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
    data_long <- data_wide[, c("YEAR", "MONTH", "DAY", "DATE", "GEAR_CAT", 
                               "MEAN_LONGITUDE", "MEAN_LATITUDE", "HAULJOIN", 
                               "GEAR_DEPTH", "GEAR_TEMPERATURE", "WEIGHT",
                               "AREA_SWEPT", ispp_code)]
    names(data_long) <- c("year", "month", "day", "date", "gear", "lon", "lat", 
                          "hauljoin", "bot_depth", "bot_temp", "catch_kg",
                          "area_swept_km2", "cpue_kg_km2")
    
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
                                     "GEAR_CAT", "MEAN_LONGITUDE", 
                                     "MEAN_LATITUDE", "HAULJOIN", 
                                     "GEAR_DEPTH", "GEAR_TEMPERATURE", 
                                     "WEIGHT", "AREA_SWEPT")],
                       cpue_kg_km2 = rowSums(sub_df))
    names(data_long) <- c("year", "month", "day", "date", "gear", "lon", "lat",
                          "hauljoin", "bot_depth", "bot_temp", "catch_kg",
                          "area_swept_km2", "cpue_kg_km2")
    
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