##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rmarkdown)
library(tidyr)
library(kableExtra)
library(FishStatsUtils)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set draft version
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draft_version <- 4
draft_dir <- paste0("manuscript/Version_", draft_version, "/")
if(!dir.exists(draft_dir)) dir.create(draft_dir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Copy the most current version of the figures into the draft version folder
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file.copy(from = "figures/",
          to = draft_dir,
          recursive = TRUE)

n_figs <- length(x = grep(x = dir("figures/"), pattern = ".jpeg"))
n_tables <- 1
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Render the main manuscript text
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

todays_date <- paste0(format(x = Sys.Date(), "%d"), 
                      format(x = Sys.Date(), "%m"), 
                      format(x = Sys.Date(), "%y"), 
                      collapse = "")

rmarkdown::render(input = "manuscript/Oyafuso_etal_Chukchi_MS.Rmd", 
                  output_format = "word_document", 
                  output_file = paste0("../", draft_dir,
                                       "Oyafuso_Chukchi Survey Designv", 
                                       draft_version, "_", todays_date,".docx"),
                  params = list(draft_dir = draft_dir, 
                                n_figs = n_figs, n_tables = n_tables))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Render Appendix A and B
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rmarkdown::render(input = "manuscript/Oyafuso_etal_Chukchi_Appendix_A.Rmd",
#                   output_format = "word_document",
#                   output_dir =  draft_dir,
#                   output_file = "Appendix_A.docx")
# 
# rmarkdown::render(input = "manuscript/Oyafuso_etal_Chukchi_Appendix_B.Rmd",
#                   output_format = "word_document",
#                   output_dir = draft_dir,
#                   output_file = "Appendix_B.docx",
#                   params = list(draft_dir = draft_dir))
