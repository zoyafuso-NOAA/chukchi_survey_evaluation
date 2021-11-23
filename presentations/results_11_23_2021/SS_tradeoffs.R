load("data/survey_opt_data/optimization_data_otter.RData")

par(mfrow = c(3, 2), mar = c(3, 3, 1, 1))
for (ispp in spp_list) {
  n_runs <- length(dir(paste0("results/otter_trawl/survey_opt/SS/", ispp, "/")))
  
  cv_n <- data.frame(run = 1:n_runs)
  
  for (irun in 1:n_runs){
    load(paste0("results/otter_trawl/survey_opt/SS/", ispp, 
                "/Run_", irun, "/result_list.RData"))
    cv_n$cv[irun] <- result_list$cvs
    cv_n$n[irun] <- result_list$n
  }
  
  plot(cv ~ n, data = cv_n, pch = 16, las = 1, ann = F)
  lines(cv ~ n, data = cv_n)
  legend("topright", legend = ispp, bty = "n", text.font = 2)
  
}

