##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Synthesize Alaska Bottom Trawl Arctic beam trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(dplyr)
library(readxl)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
main_haul <- as.data.frame(readxl::read_xlsx(
  "data/fish_data/2017_2019_Beam/Hauls2017_2019.xlsx"))
main_haul <- subset(main_haul, 
                      select = c("Year" , "EVENT_EVENT_NUMBER",
                                 "BOTTOM_DEPTH_m", "AREA_SWEPT_KM_2",
                                 "TIMESTAMP_ON_BOTTOM_GMT",
                                 "LAT_DD_ON_BOTTOM", "LON_DD_ON_BOTTOM",
                                 "Surface_temp_C", "Surface_sal", 
                                 "bottom_temp_C", "bottom_sal") )

main_catch <- as.data.frame(readxl::read_xlsx(
  "data/fish_data/2017_2019_Beam/Genetics_corrected_Catch.xlsx"))

species_codes <- read.csv("data/fish_data/AK_BTS_OtterAndBeam/species.csv")
taxa <- read.csv("data/fish_data/AK_BTS_OtterAndBeam/taxonomy.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create  result directory ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_dir <- "data/fish_data/2017_2019_Beam/data_long_by_taxa/"
if (!dir.exists(output_dir)) dir.create(output_dir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Subset species ----
##   Subset catch data to selected fish species 
##   Sum over subsampling?
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main_catch <- subset(x = main_catch,
                       select = c(CLAMS_EVENT_NUMBER, 
                                  TOTAL_WEIGHT_IN_HAUL,
                                  SPECIES_CODE))
names(main_catch) <- c("station_id", "catch_kg", "species_code")

main_catch <- aggregate(catch_kg ~ species_code + station_id, 
                          data = main_catch,
                          FUN = sum)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Match haul data ----
##   Match haul data information with catch data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
idx <- match(main_catch$station_id,
             main_haul$EVENT_EVENT_NUMBER)

main_catch[, c("year", "station_id", "area_swept_km2", "date",
               "lat", "lon", "bot_temp", "bot_depth")] <-
  main_haul[idx, c("Year", "EVENT_EVENT_NUMBER", "AREA_SWEPT_KM_2",
                   "TIMESTAMP_ON_BOTTOM_GMT",
                   "LAT_DD_ON_BOTTOM", "LON_DD_ON_BOTTOM", 
                   "bottom_temp_C", "BOTTOM_DEPTH_m")]
main_catch$cpue_kg_km2 <- main_catch$catch_kg / main_catch$area_swept_km2

main_catch$date <- as.Date(x = main_catch$date, format = "%d-%b-%y")
main_catch$month <- lubridate::month(x = main_catch$date)
main_catch$day <- lubridate::day(x = main_catch$date)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Widen data ----
##   Create a wide df with zeros filled for stations where species not observed
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_wide <- tidyr::spread(data = subset(main_catch,
                                         select = -catch_kg),
                           key = "species_code", 
                           value = "cpue_kg_km2",
                           fill = 0)
data_wide <- data_wide[!is.na(data_wide$year), ]
data_wide$gear <- "beam"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create datasets ----
##   Create VAST data input for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_info <- data.frame()

solo_spp <- subset(x = species_codes, 
                   subset = COMMON_NAME %in% 
                     c("Alaska plaice", "Arctic cod",  "Pacific cod", 
                       "walleye pollock", "snow crab", "yellowfin sole", 
                       "Bering flounder", "saffron cod", 
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
    data_long <- cbind(data_wide[, c("station_id", "year", "month",
                                     "day", "date", "area_swept_km2", "gear",
                                     "lat", "lon", "bot_temp", "bot_depth")],
                       cpue_kg_km2 = data_wide[, paste(ispp_code)])

    
    ## Remove species column from data_wide
    data_wide <- data_wide[!names(data_wide) %in% paste(ispp_code)]
    
    ## Save data
    write.csv(x = data_long, file = paste0(output_dir, ispp, ".csv"), 
              row.names = F) 
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Aggregate species ----
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
    
    data_long <- cbind(data_wide[, c("station_id", "year", "month",
                                     "day", "date", "area_swept_km2", "gear",
                                     "lat", "lon", "bot_temp", "bot_depth")],
                       cpue_kg_km2 = rowSums(sub_df))
    
    ## Remove species column from data_wide
    data_wide <- data_wide[!names(data_wide) %in% paste(ispp_code)]
    
    ## Save data
    write.csv(x = data_long, file = paste0(output_dir, ispp, ".csv"),
              row.names = F) 
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save ----
##   Create VAST data input for each species group
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
write.csv(x = species_info, 
          file = paste0(dirname(output_dir), "/species_info.csv"),
          row.names = F)