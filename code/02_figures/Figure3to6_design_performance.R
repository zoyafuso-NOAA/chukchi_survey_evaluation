##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Simulation Performance 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot, True CV, RRMSE of CV, and Bias, for SRS, STRS, and SYS 
##                    designs for each species.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Plot True CV and RRMSE of CV ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (igear in c("otter", "beam")) {
  
  ## Load data ----
  load(file = paste0("results/chukchi_", igear, 
                     "/random_survey_sim_results.RData"))
  
  ## Calculate sample sizes for systematic grid
  sys_n <- aggregate(n ~ res, data = sys_settings, FUN = mean)
  sys_n <- sys_n[order(sys_n$res, decreasing = TRUE), ]

  spp_name <- dimnames(x = true_cv_srs)[[2]]
  
  for (iplot in c("performance", "bias")) {
    fig.number <- switch(paste0(igear, "_", iplot),
                         "otter_performance" = 3,
                         "otter_bias" = 4,
                         "beam_performance" = 5,
                         "beam_bias" = 6)
    
    file_name <- paste0("figures/Figure", fig.number, "_", 
                        igear, "_", iplot, ".png")
    metrics <- 
      list(performance = c("True CV" = "true_cv", 
                           "RRMSE of CV" = "rrmse_cv"),
           bias = c("% Bias of Index" = "rb_index",
                    "% Bias of CV Relative to True CV" = "rb_cv"))[[iplot]]
    
    #   Open device
    png(filename = file_name, res = 500, family = "serif", units = "in", 
        width = 19.0 / 2.54, 
        height = c("otter" = 16 / 2.54, "beam" = 16 / 2.54)[igear] )
    
    ##   Plot layout 
    par(mfrow = list("otter" = c(4, 8), 
                     "beam" = c(4, 8))[[igear]], 
        oma = c(2, 3, 0, 3))
    
    for (ispp in 1:length(spp_name)) { ## Loop over species
      for (imetric in 1:length(metrics)) {
        
        temp_metric <- metrics[imetric]
        par(mar = switch(imetric, "1" = c(2, 2, 2, 0), 
                         "2" = c(2, 0, 2, 2)))
        
        if (iplot == "bias") {# Calculate mean for ispp across sample sizes
          srs_metric <- apply(X = get(paste0(temp_metric[1], "_srs")), 
                              MARGIN = list(c(1, 2), c(2, 3))[[imetric]], 
                              FUN = "mean", 
                              na.rm = TRUE)[, ispp]
          strs_metric <- apply(X = get(paste0(temp_metric[1], "_ms_strs")), 
                               MARGIN = list(c(1, 2), c(2, 3))[[imetric]], 
                               FUN = "mean", 
                               na.rm = TRUE)[, ispp]
          sys_metric <- apply(X = get(paste0(temp_metric[1], "_sys")), 
                              MARGIN = list(c(1, 2), c(2, 3))[[imetric]], 
                              FUN = "mean", 
                              na.rm = TRUE)[, ispp]
        }
        
        if (iplot == "performance") {
          srs_metric <- get(paste0(temp_metric[1], "_srs"))[, ispp]
          strs_metric <- get(paste0(temp_metric[1], "_ms_strs"))[, ispp]
          sys_metric <- get(paste0(temp_metric[1], "_sys"))[, ispp]
        }
        
        ylim_ <- switch(iplot, 
                        "performance" = {
                          bound = max(srs_metric, strs_metric, sys_metric, 
                                      na.rm = TRUE)
                          c(0, 1.1 * bound)},
                        "bias" = {
                          bound <- max(5, abs(c(srs_metric, 
                                                 strs_metric, 
                                                 sys_metric)), na.rm = TRUE)
                          c(-1.2 * bound, 1.2 * bound)
                        }
        )
        
        ## Base Plot
        plot(1, type = "n", axes = F, ann = F, 
             ylim = ylim_, 
             xlim = c(45, 180))
        
        ## add simple and stratified random metrics
        matpoints(x = target_n, 
                  y = cbind(srs_metric, strs_metric), 
                  pch = 16, col = c("black", "orange"))
        matlines(x = target_n, 
                 y = cbind(srs_metric, strs_metric), 
                 lty = 1, col = c("black", "orange"))
        points(x = sys_n$n, y = sys_metric, 
               pch = 16, col = "blue")
        lines(x = sys_n$n, y = sys_metric, 
              pch = 16, col = "blue")
        
        axis(side = 1, cex.axis = 0.8, at = c(50, 100, 150))
        axis(side = ifelse(test = imetric == 1, yes = 2, no = 4), 
             las = 1, cex.axis = 1, tick = F,
             line = ifelse(test = imetric == 1, yes = -0.5, no = -0.5),  
             labels = pretty(ylim_, n = 3),
             at = pretty(ylim_, n = 3))
        axis(side = ifelse(test = imetric == 1, yes = 2, no = 4), 
             las = 1, cex.axis = 1, tick = T, tcl = -0.25,
             labels = NA,
             at = pretty(ylim_, n = 3))
        box()
        if (iplot == "bias") abline(h = 0)
        
        ## Plot species subtitle
        if (imetric == 1)
          text(x = 200,
               y = ylim_[2] * c("performance" = 1.15, "bias" = 1.3)[iplot],
               labels = spp_name[ispp],
               font = 2,
               xpd = NA)
        
      }
    }
    
    ##  Plot legend
    par(mar = c(2, 2, 2, 0))
    plot(1, type = "n", axes = F, ann = F,
         ylim = c(0, 100), xlim = c(0, 100))
    legend(x = 5, y = 90, bty = "n", cex = 1.5,
           legend = c("SRS", "MS-STRS", "SYS"),
           title = paste(igear, "trawl"),
           col = c("black", "orange", "blue"),
           lty = 1, pch = 16, xpd = NA)
    
    ## Axis Labels 
    mtext(side = 2, text = names(metrics)[1], 
          outer = TRUE, line = 1, font = 2)
    mtext(side = 4, text = names(metrics)[2], 
          outer = TRUE, line = 1, font = 2)
    mtext(side = 1, text = "Sample Size", 
          outer = TRUE, line = 1, font = 2)
    dev.off()
  }
}
