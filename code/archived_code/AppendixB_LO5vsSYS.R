##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Simulation Performance: SYS estimator comparison
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare the RRMSE of CV and bias of CV of the conventional
##                      systematic variance estimator versus the localize
##                      five-neighbor (LO5) estimator
##
## Notes:         
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

for (igear in c("beam", "otter")) { ## loop over gear -- start
  
  ## Load survey simulation results for a given igear
  load(paste0("results/chukchi_", igear, "/random_survey_sim_results.RData"))
  
  ## Get species list for a given igear
  spp_list <- dimnames(x = rrmse_cv_sys)[[2]]
  
  ## Plot file name
  file_name <- paste0("figures/AppendixB_", igear, ".png")
  
  ## Open png file
  png(filename = file_name, res = 500, family = "serif", units = "in", 
      width = 19.0 / 2.54, 
      height = c("otter" = 16 / 2.54, "beam" = 16 / 2.54)[igear] )
  
  ##   Plot layout 
  par(mfrow = list("otter" = c(4, 8), 
                   "beam" = c(4, 8))[[igear]], 
      oma = c(2, 3, 0, 3))
  
  for (ispp in spp_list) { ## loop over species -- start
    for (imetric in c("rrmse", "rb_cv")) { ## loop over metric -- start
      
      ## Modify plot margins based on the metric. This makes it so that the 
      ## two plots touch each other side by side.
      par(mar = switch(imetric, "rrmse" = c(2, 2, 2, 0), 
                       "rb_cv" = c(2, 0, 2, 2)))
      
      ## Extract the imetric for a given ispp
      plot_this <- 
        switch(imetric,
               "rrmse" = cbind(rrmse_cv_sys[, ispp],
                               rrmse_cv_sys_LO5[, ispp]),
               "rb_cv" = cbind(apply(rb_cv_sys, 
                                     MARGIN = 2:3, 
                                     FUN = mean, 
                                     na.rm = TRUE)[, ispp],
                               apply(rb_cv_sys_LO5, 
                                     MARGIN = 2:3, 
                                     FUN = mean, 
                                     na.rm = TRUE)[, ispp]))
      
      ## Calculate the ylim of the plot based on the imetric. Bias plots
      ## are symetrical in the y-direction.
      ylim_ <- switch(imetric, 
                      "rrmse" = {
                        bound = max(plot_this, na.rm = TRUE)
                        c(0, 1.1 * bound)},
                      "rb_cv" = {
                        bound <- max(5, abs(plot_this), na.rm = TRUE)
                        c(-1.2 * bound, 1.2 * bound)
                      }
      )
      
      ## Plot
      matplot(x = sys_settings$n,
              y = plot_this,
              xlim = c(45, 180), ylim = ylim_,
              axes = F,
              pch = 16, col = c("blue", "orange"), las = 1 )
      matlines(x = sys_settings$n,
               y = plot_this,
               col = c("blue", "orange"), lty = 1)
      
      ## add axis and tick labels
      axis(side = 1, cex.axis = 0.65, at = c(50, 100, 150))
      axis(side = c("rrmse" = 2, "rb_cv" = 4)[imetric], 
           las = 1, cex.axis = 0.7, 
           labels = pretty(ylim_, n = 3),
           at = pretty(ylim_, n = 3))
      box()
      
      if (imetric == "rb_cv") abline(h = 0, lty = 'dotted')
      
      ## Plot species subtitle
      if (imetric == "rrmse")
        text(x = 200,
             y = ylim_[2] *  1.15,
             labels = ispp,
             font = 2,
             xpd = NA)
    }  ## loop over metric -- end
  } ## loop over species -- end
  
  ##  Plot legend
  par(mar = c(2, 2, 2, 0))
  plot(1, type = "n", axes = F, ann = F,
       ylim = c(0, 100), xlim = c(0, 100))
  legend(x = 5, y = 90, bty = "n", cex = 1.5,
         legend = c("SYS", "SYS-LO5"),
         title = paste(igear, "trawl"),
         col = c("blue", "orange"),
         lty = 1, pch = 16, xpd = NA)
  
  ## Axis Labels 
  mtext(side = 2, text = "RRMSE of CV", 
        outer = TRUE, line = 1, font = 2)
  mtext(side = 4, text = "% Bias of CV Relative to True CV", 
        outer = TRUE, line = 1, font = 2)
  mtext(side = 1, text = "Sample Size", 
        outer = TRUE, line = 1, font = 2)
  
  ## Close plot device
  dev.off()
} ## loop over gear -- end
