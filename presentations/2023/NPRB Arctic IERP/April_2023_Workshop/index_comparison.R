##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Preliminary Plots for
##   Arctic IERP April 2023 Workshop
##   Abundance indices for NBS / Chukchi taxa (otter trawl only?)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Authorize Google Drive
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(googledrive)
googledrive::drive_deauth()
googledrive::drive_auth()
1

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants: googledrive ids point to folder in the AIERP Synthesis 
##   Living Documents / Preliminary Analyses
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbs_otter_dir <- googledrive::as_id("1eERDrW32SGU9hABm7bRC7y5gnUAUafeg")
chukchi_otter_dir <- googledrive::as_id("1J1UW1Q_UkO-tVaATpm_nWYdXfw6Rijvw")

nbs_otter_ls <- googledrive::drive_ls(path = nbs_otter_dir)
chukchi_otter_ls <- googledrive::drive_ls(path = chukchi_otter_dir)

nbs_spp <- sort(nbs_otter_ls$name)
chukchi_spp <- sort(chukchi_otter_ls$name)

nbs_otter_years <- c(1985, 1988, 1991, 2010, 2017:2019, 2021:2022)
chukchi_otter_years <- c(1990, 2012)
chukchi_beam_years <- c(2012, 2017, 2019)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Extract abundance indices for each taxon / region
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

index <- data.frame() ## main df 

for (ispp in nbs_spp) {
  which_surveys <- c("nbs_otter", "chukchi_otter")[c(ispp %in% nbs_spp, 
                                                     ispp %in% chukchi_spp)] 
  
  for (isurvey in which_surveys) {
    
    ## Find folder where Index.csv lives for a given taxon / region
    temp_dir_df <- get(x = paste0(isurvey, "_ls"))
    temp_spp <- temp_dir_df$id[temp_dir_df$name == ispp]
    temp_years <- get(x = paste0(isurvey, "_years"))
    
    diagnostics_id <- googledrive::drive_ls(path = temp_spp,
                                            pattern = "diagnostics")$id
    
    if (length(x = diagnostics_id) == 0) next
    
    ## Find file that holds Index.csv for a given taxon/region
    index_id <- googledrive::drive_ls(path = diagnostics_id, 
                                      pattern = "Index.csv")$id
    index_string <- googledrive::drive_read_string(file = index_id, 
                                                   encoding = "UTF-8")
    
    ## Import Index.csv
    temp_index <- read.csv(text = index_string)
    temp_index <- subset(x = temp_index, 
                         subset = Time %in% temp_years, 
                         select = c(Time, Estimate, Std..Error.for.Estimate, 
                                    Std..Error.for.ln.Estimate.))
    
    ## Append to main df
    index <- rbind(index, 
                   cbind("Survey" = isurvey, 
                         "Taxon" = ispp,
                         temp_index))
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

index$lower_sd <- index$Estimate - index$Std..Error.for.Estimate
index$upper_sd <- index$Estimate + index$Std..Error.for.Estimate

index$Estimate <- index$Estimate / 1e6
index$lower_sd <- index$lower_sd / 1e6
index$upper_sd <- index$upper_sd / 1e6

for (ispp in nbs_spp[!nbs_spp %in% "urchins"]) {
  which_surveys <- c("nbs_otter", "chukchi_otter")[c(ispp %in% nbs_spp, 
                                                     ispp %in% chukchi_spp)] 
  
  temp_df <- subset(x = index, subset = Taxon == ispp)
  
  plot(1, type = "n", 
       xlim = c(1985, 2022),
       ylim = c(0, max(temp_df$upper_sd * 1.01)),
       yaxs = "i",
       las = 1, xlab = "", ylab = "")
  
  for (isurvey in which_surveys) {
    points(Estimate ~ Time, data = temp_df, subset = Survey == isurvey, 
           pch = 20, 
           col = c("nbs_otter" = "black", "chukchi_otter" = "blue")[isurvey])
    lines(Estimate ~ Time, data = temp_df, subset = Survey == isurvey, 
           col = c("nbs_otter" = "black", "chukchi_otter" = "blue")[isurvey])
    segments(x0 = temp_df$Time[temp_df$Survey == isurvey],
             x1 = temp_df$Time[temp_df$Survey == isurvey],
             y0 = temp_df$lower_sd[temp_df$Survey == isurvey],
             y1 = temp_df$upper_sd[temp_df$Survey == isurvey], 
             col = c("nbs_otter" = "black", "chukchi_otter" = "blue")[isurvey])
  }
  
  mtext(side = 3, text = ispp)
}
