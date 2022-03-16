###############################################################################
## Project:       Simulation Performance 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot the True CV and RRMSE of CV for SRS, STRS, and SYS 
##                    designs for each species.
###############################################################################
rm(list = ls())

##################################################
####   Choose a gear
##################################################  
igear = "otter"

##################################################
####   Load data
##################################################  
load(paste0("results/chukchi_", igear, "/survey_sim_results.RData"), 
     verbose = TRUE)

settings <- list(
  otter = list(
    spp_idx = list(fish = c(1:3, 5, 7, 9, 10, 12, 13, 15, 19, 20),
                   inverts = c(4, 6, 8, 11, 14, 16:18)),
    mfrow = list(fish = c(4, 6),
                 inverts = c(3, 6)),
    plot_height_cm = c(fish = 20, inverts = 15) 
  )
)


for (taxon in c("fish", "inverts")) {
  
  plot_height_in <- settings[[igear]][["plot_height_cm"]][[taxon]] / 2.54
  plot_mfrow <- settings[[igear]][["mfrow"]][[taxon]]
  idx <- settings[[igear]][["spp_idx"]][[taxon]]
  spp_name <- dimnames(true_cv_srs)[[2]][idx]
  
  ##################################################
  ####   Open device
  ##################################################  
  png(filename = paste0("figures/FigureX_", igear, 
                        "_", taxon, "_performance.png"), 
      width = 19.0 / 2.54, 
      height = plot_height_in,
      res = 500, units = "in")

  ##################################################
  ####   Plot layout
  ##################################################  
  par(mfrow = plot_mfrow,
      oma = c(3, 3, 0, 3))
  
  for (ispp in 1:length(idx)) { ## Loop over species -- start
    
    ##################################################
    ####   Plot True CV
    ##################################################  
    par(mar = c(2.5, 2.5, 2.5, 0))
    ylim = max(true_cv_srs[, idx[ispp] ],
               true_cv_sys[, idx[ispp]],
               true_cv_ms_strs[, idx[ispp]]) * 1.1
    
    plot(1, type = "n", axes = F, ann = F, 
         ylim = c(0, ylim), xlim = c(35, 180))
    
    points(x = target_n, y = true_cv_srs[, idx[ispp]], pch = 16)
    lines(x = target_n, y = true_cv_srs[, idx[ispp]], pch = 16)
    points(x = sys_settings$n, y = true_cv_sys[, idx[ispp]], 
           pch = 16, col = "blue")
    lines(x = sys_settings$n, y = true_cv_sys[, idx[ispp]], 
          pch = 16, col = "blue")
    points(x = target_n, y = true_cv_ms_strs[, idx[ispp]], 
           pch = 16, col = "orange")
    lines(x = target_n, y = true_cv_ms_strs[, idx[ispp]], 
          pch = 16, col = "orange")

    axis(side = 1, cex.axis = 0.75)
    axis(side = 2, las = 1, cex.axis = 0.75)
    box()

    ##################################################
    ####   Design Legend
    ##################################################  
    if (ispp %in% c(1, length(idx))) {
      legend("bottom", bty = "n", cex = 0.75,
             legend = c("SRS", "MS-STRS", "SYS"),
             col = c("black", "orange", "blue"),
             lty = 1, pch = 16)
    }
    
    ##################################################
    ####   Plot Species Subtitle
    ##################################################
    text(x = 110,
         y = ylim * 1.15,
         labels = spp_name[ispp],
         font = 2,
         xpd = NA)

    ##################################################
    ####   Plot RRMSE of CV
    ##################################################
    par(mar = c(2.5, 0, 2.5, 2.5))
    ylim = max(rrmse_cv_srs[, idx[ispp]],
               rrmse_cv_sys[, idx[ispp]],
               rrmse_cv_ms_strs[, idx[ispp]]) * 1.1
    
    plot(1, type = "n", axes = F, ann = F,
         ylim = c(0, ylim), xlim = c(35, 180))
    
    points(x = target_n, y = rrmse_cv_srs[, idx[ispp]], pch = 16)
    lines(x = target_n, y = rrmse_cv_srs[, idx[ispp]], pch = 16)
    points(x = sys_settings$n, y = rrmse_cv_sys[, idx[ispp]], 
           pch = 16, col = "blue")
    lines(x = sys_settings$n, y = rrmse_cv_sys[, idx[ispp]], 
          pch = 16, col = "blue")
    points(x = target_n, y = rrmse_cv_ms_strs[, idx[ispp]], 
           pch = 16, col = "orange")
    lines(x = target_n, y = rrmse_cv_ms_strs[, idx[ispp]], 
          pch = 16, col = "orange")

    axis(side = 1, cex.axis = 0.75)
    axis(side = 4, las = 1, cex.axis = 0.75)
    box()
    
  }
  
  ## Axis Labels
  mtext(side = 2, text = "True CV", outer = TRUE, line = 1, font = 2)
  mtext(side = 4, text = "RRMSE of CV", outer = TRUE, line = 1, font = 2)
  mtext(side = 1, text = "Sample Size", outer = TRUE, line = 1, font = 2)
  dev.off()
} ## Loop over species -- end
