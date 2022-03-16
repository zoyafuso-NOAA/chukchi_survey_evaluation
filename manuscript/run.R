library(rmarkdown)

draft_version <- 0
draft_dir <- paste0("manuscript/Version_", draft_version, "/")
if(!dir.exists(draft_dir)) dir.create(draft_dir)

rmarkdown::render(input = "manuscript/Oyafuso_etal_Chukchi_survey_testing.Rmd", 
                  output_format = "word_document", 
                  output_file = paste0("../", draft_dir,
                                      "Oyafuso_etal_Chukchi_survey_testing.docx"),
                  params = list(draft_dir = draft_dir))

file.copy(from = "figures/",
          to = draft_dir,
          recursive = TRUE)

