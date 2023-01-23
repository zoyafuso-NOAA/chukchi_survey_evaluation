###############################################################################
## Project:       Ckean up single-species 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################
## Import constants
##################################
load("data/survey_opt_data/optimization_data.RData")

##################################
## Result objects
##################################
for (igear in c("beam", "otter")) {
  n_spp <- get(paste0("n_spp_", igear))
  spp_list <- get(paste0("spp_list_", igear))
  
  ss_solutions <- matrix(nrow = n_cells, ncol = n_spp,
                         dimnames = list(NULL, spp_list))
  
  ss_strs_cvs <- vector(length = n_spp)
  target_n <- 100
  
  ss_allocations <- ss_population <-  
    matrix(data = 0, nrow = 5, ncol = n_spp,
           dimnames = list(paste0("Str_", 1:5), spp_list))
  
  ##################################
  ## 
  ##################################
  
  for (ispp in 1:n_spp) {
    n_runs <- length(dir(paste0("results/",
                                ifelse(test = igear == "otter",
                                       yes = "otter_trawl",
                                       no = "beam_trawl_2012_2019"),
                                "/survey_opt/SS/", 
                                spp_list[ispp], "/")))
    cv_n <- data.frame(run = 1:n_runs)
    
    for (irun in 1:n_runs){
      load(paste0("results/",
                  ifelse(test = igear == "otter",
                         yes = "otter_trawl",
                         no = "beam_trawl_2012_2019"),
                  "/survey_opt/SS/", spp_list[ispp], 
                  "/Run_", irun, "/result_list.RData"))
      cv_n$cv[irun] <- result_list$cvs
      cv_n$n[irun] <- result_list$n
    }
    
    idx <- which.min(abs(cv_n$n - target_n))
    cv_n[idx, ]
    
    load(paste0("results/",
                ifelse(test = igear == "otter",
                       yes = "otter_trawl",
                       no = "beam_trawl_2012_2019"),
                "/survey_opt/SS/", spp_list[ispp], 
                "/Run_", idx, "/result_list.RData"))
    
    error_df <- data.frame("DOM" = "DOM1",
                           "CV1" = result_list$cvs,
                           "domainvalue"  = 1)
    
    temp_stratif <- result_list$solution$aggr_strata
    temp_stratif$DOM1 <- 1
    
    temp_bethel <- SamplingStrata::bethel(
      errors = error_df,
      stratif = temp_stratif, 
      realAllocation = T, 
      printa = T)
    temp_n <- sum(ceiling(temp_bethel))
    
    while (temp_n != target_n){
      over_under <- temp_n > target_n
      CV_adj <- ifelse(over_under == TRUE, 
                       yes = 1.001,
                       no = 0.999)
      
      error_df$CV1 <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * CV_adj
      temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                            errors = error_df, 
                                            printa = TRUE)
      
      temp_n <- sum(as.numeric(temp_bethel))
      
      print(paste0("n = ", temp_n, ", CV = ", 
                   as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])) )
    }
    
    ss_allocations[1:length(temp_bethel), ispp] <- round(as.numeric(temp_bethel))
    ss_strs_cvs[ispp] <- as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
    ss_solutions[, ispp] <- result_list$sol_by_cell
    ss_population[1:length(temp_bethel), ispp] <- result_list$sum_stats$Population
    
    for (ivar in c("ss_allocations", "ss_strs_cvs", 
                   "ss_solutions", "ss_population")){
      assign(x = paste0(ivar, "_", igear), value = get(ivar))
    }
    
    ##################################
    ## Save
    ##################################
    save(list = paste0(c("ss_allocations", "ss_strs_cvs", 
                         "ss_solutions", "ss_population"), "_", igear),
         file = paste0("results/",
                       ifelse(test = igear == "otter",
                              yes = "otter_trawl",
                              no = "beam_trawl_2012_2019"),
                       "/survey_opt/SS/ss_knit_results_", 
                       igear, ".RData"))
  }
  
  
}
