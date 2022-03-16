library(tools)

img_dir <- "//akc0ss-n086/AKC_PubliC/Dropbox/RACE Survey App Legacy/RACESurveyApp_2017-3/My Web Sites/Fish ID_files/"

data <- read.csv("data/taxon_images/taxon_images.csv")

for (irow in 9) {
  if(data$file_name[irow] != "") {
    file_type <- file_ext(paste0(img_dir, data$file_name[irow]))
    file.copy(from = paste0(img_dir, data$file_name[irow]),
              to = paste0("data/taxon_images/", 
                          data$common_name[irow],
                          ".", file_type))
  }  
}
