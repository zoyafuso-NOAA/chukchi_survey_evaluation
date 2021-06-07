###############################################################################
## Project:       Synthesize Alaska Bottom Trawl Arctic beam trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries
##################################################
library(dplyr)
library(readxl)

##################################################
####  Import Data
##################################################
master_haul <- as.data.frame(readxl::read_xlsx(
  "data/fish_data/2017_2019_Beam/Hauls2017_2019.xlsx"))
master_haul <- subset(master_haul, 
                      select = c("Year" , "EVENT_EVENT_NUMBER",
                                 "BOTTOM_DEPTH_m", "AREA_SWEPT_KM_2",
                                 "LAT_DD_ON_BOTTOM", "LON_DD_ON_BOTTOM",
                                 "Surface_temp_C", "Surface_sal", 
                                 "bottom_temp_C", "bottom_sal") )

master_catch <- as.data.frame(readxl::read_xlsx(
  "data/fish_data/2017_2019_Beam/Genetics_corrected_Catch.xlsx"))

##################################################
####  Subset catch data to selected fish species 
####  Sum over subsampling?
##################################################
master_catch <- subset(x = master_catch,
                       select = c(CLAMS_EVENT_NUMBER, 
                                  TOTAL_WEIGHT_IN_HAUL,
                                  COMMON_NAME))
names(master_catch) <- c("station_id", "catch_kg", "common_name")

master_catch <- aggregate(catch_kg ~ common_name + station_id, 
                          data = master_catch,
                          FUN = sum)

##################################################
####  Match haul data information with catch data
##################################################
idx <- match(master_catch$station_id,
             master_haul$EVENT_EVENT_NUMBER)

master_catch[, c("year", "station_id", "area_swept_km2", "lat", "lon")] <-
  master_haul[idx, c("Year", "EVENT_EVENT_NUMBER", "AREA_SWEPT_KM_2",
                     "LAT_DD_ON_BOTTOM", "LON_DD_ON_BOTTOM")]

##################################################
#### Create a wide df with zeros filled for stations where species not observed
##################################################
data_wide <- tidyr::spread(data = master_catch,
                           key = "common_name", 
                           value = "catch_kg",
                           fill = 0)
data_wide <- data_wide[!is.na(data_wide$year), ]

##################################################
#### "lengthen" wide dataset to long-form to match VAST data input 
##################################################
data_long <- reshape::melt(data_wide, 
                           measure.vars = unique(master_catch$common_name),
                           variable_name = "common_name")

##################################################
#### Subset data to a set of species of interest
#### Change common names to be consistent with RACEBACE
##################################################
spp_list <- c("Arctic cod", "Bering flounder", "Saffron cod", "Snow crab",
              "Pacific cod", 
              "Walleye pollock", "Yellowfin sole", "Alaska plaice",
              "Starry flounder", "Arctic staghorn sculpin", "Pacific herring",
              "Slender eelblenny", "Blue king crab", "Hairy hermit crab",
              "Shorthorn (Warty) sculpin", "Sculptured shrimp", 
              "Notched brittlestar", "Purple Orange sea star", 
              "Fuzzy hermit crab", "Circumboreal toad crab", "Humpy shrimp",
              "Helmet crab", "Arctic argid", "Kuro argid", "Green sea urchin",
              "Northern nutclam", "Common mud star", "Greenland cockle",
              "Basketstar")

# c("rainbow smelt", "circumpolar eualid", "shortscale eualid", 
# "smooth nutclam", "brownscaled sea cucumber")

data_long <- subset(x = data_long, 
                    subset = common_name %in% spp_list)
data_long$common_name <- as.character(data_long$common_name)

data_long$common_name <- unlist(sapply(
  X = data_long$common_name,
  FUN = function(x) switch(
    x,
    "Arctic cod" = "Arctic cod", 
    "Bering flounder" = "Bering flounder", 
    "Starry flounder" = "starry flounder",
    "Saffron cod" = "saffron cod", 
    "Walleye pollock" = "walleye pollock",
    "Pacific cod" = "Pacific cod",
    "Yellowfin sole" = "yellowfin sole",
    "Alaska plaice" = "Alaska plaice", 
    "Arctic staghorn sculpin" = "Arctic staghorn sculpin", 
    "Pacific herring" = "Pacific herring"  , 
    "Slender eelblenny" = "slender eelblenny", 
    "Shorthorn (Warty) sculpin" = "shorthorn (=warty) sculpin", 
    
    "Snow crab" = "snow crab",
    "Blue king crab" = "blue king crab", 
    "Fuzzy hermit crab" = "fuzzy hermit crab", 
    "Circumboreal toad crab" = "circumboreal toad crab", 
    "Helmet crab" = "helmet crab",
    "Hairy hermit crab" = "hairy hermit crab",
    
    "Notched brittlestar" = "notched brittlestar", 
    "Purple Orange sea star" = "purple-orange sea star", 
    "Basketstar" = "basketstar",
    "Common mud star" = "common mud star", 
    "Green sea urchin" = "green sea urchin", 
    
    "Humpy shrimp" = "humpy shrimp",
    "Sculptured shrimp" = "sculptured shrimp", 
    "Arctic argid" = "Arctic argid", 
    "Kuro argid" = "kuro argid", 
    
    "Northern nutclam" = "northern nutclam", 
    "Greenland cockle" = "Greenland cockle"
  )))

##################################################
####  Prepare the dataframe for catch-rate data in the VAST format
##################################################
ierl_data_long <- subset(x = data_long, 
                         select = -station_id)
names(ierl_data_long)[6] <- "catch_kg"
ierl_data_long$gear <- "beam"

##################################################
####  Save
##################################################
write.csv(x = ierl_data_long, 
          file = "data/fish_data/2017_2019_Beam/ierl_data_processed.csv", 
          row.names = FALSE)
