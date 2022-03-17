###############################################################################
## Project:       Figure X: VAST output
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   predicted density, index, qq plots for each species/gear
###############################################################################
rm(list = ls())

##################################################
#### Import relevant libraries   
##################################################
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)
library(gap)
library(plotrix)
library(VAST)
library(jpeg)

##################################################
#### Import Data and Constants
##################################################
load("data/survey_opt_data/optimization_data.RData")
AK_land <- rgdal::readOGR(dsn = "data/spatial_data/land_shapefiles/AKland.shp")
cropped_extent <- extent(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- raster::crop(x = AK_land, y = cropped_extent)
good_species <- read.csv(file = "results/good_species.csv")

##################################################
####   Result directory
##################################################  
if(!dir.exists("figures/Appendix_A/"))
  dir.create("figures/Appendix_A/")
##################################################
#### Synthesize density predictions for each gear and species
##################################################
for (igear in c("otter", "beam")) {
  
  ## Years that gear was used
  years <- list("otter" = data.frame(idx = c(1, 23),
                                     year = c(1990, 2012)),
                "beam" = data.frame(idx = c(1, 6, 8),
                                    year = c(2012, 2017, 2019)))[[igear]]
  
  ## Species that gear sampled
  species <- subset(x = read.csv("results/good_species.csv"),
                    subset = gear == igear)$taxon
  
  ## Temporary result object
  ms_dens <- array(data = NA,
                   dim = c(n_cells, length(species), nrow(years)),
                   dimnames = list(NULL, species, NULL))
  
  for (ispp in species) {
    ## Load VAST results 
    load(paste0("results/chukchi_", igear, 
                "/vast_fits/", ispp, "/fit.RData"))
    
    ## Append to temporary result object
    ms_dens[, ispp, ] <- fit$Report$D_gct[, , years$idx]
  }
  
  ## attach gear type with ms_dens variable name
  assign(x = paste0("ms_dens_", igear), value = ms_dens)
}

##################################################
#### Plot
##################################################
col_otter <- brewer.pal(n = 9, name = "Blues")
col_beam <- brewer.pal(n = 9, name = "Greens")

for (ispp_name in unique(good_species$taxon)[]){
  
  available_gears <- subset(x = good_species,
                            subset = taxon == ispp_name)$gear
  
  n_gears <- length(available_gears)
  
  png(filename = paste0("figures/Appendix_A/VAST_output_", ispp_name, ".png"),
      width = 190, units = "mm", res = 500, family = "serif",
      height = switch(paste0(available_gears, collapse = ""),
                      "beam" =  40, "otter" = 45,
                      "beamotter" = 75))
  ## Setup layout of 
  par(mar = c(0, 0, 0, 0))
  
  plot_layout_mat <- switch(paste0(available_gears, collapse = ""),
                            "beamotter" = matrix(data = c(2, 3, 6, 1, 4, 5,
                                                          12, 7, 8, 9, 10, 11),
                                                 byrow = TRUE, nrow = 2),
                            "beam" = matrix(data = c(1, 2, 3, 4, 7, 5, 6),
                                            byrow = TRUE, nrow = 1),
                            "otter" = matrix(data = c(1, 2, 3, 6, 4, 5),
                                            byrow = TRUE, nrow = 1))
  layout_widths <- switch(paste0(available_gears, collapse = ""),
                          "beamotter" = c(1, 1, 1, 1, 1.5, 1.5),
                          "beam" = c(0.5, 1, 1, 1, 0.5, 1, 1),
                          "otter" = c(1, 1, 1, 0.5, 1, 1))
  
  layout(mat = plot_layout_mat,
         widths = layout_widths)
  
  
  par(mar = c(0, 0, 0, 0))
  
  ## Name 
  plot(1, axes = F, ann = F, pch = 16, type = "n",
       xlim = c(0, 1), ylim = c(0, 1))
  mtext(side = 3, 
        text = gsub(x = ispp_name,
                    pattern = " ",
                    replacement = "\n"), 
        line = -3)
  
  if(file.exists(paste0("data/taxon_images/", ispp_name, ".jpg"))) {
    fish_img <- jpeg::readJPEG(source = paste0("data/taxon_images/", 
                                               ispp_name, ".jpg"))
    img_orientation <- ifelse(test = dim(fish_img)[1] > dim(fish_img)[2],
                              yes = "portrait", no = "landscape")
    rasterImage(image = fish_img, 
                xleft = ifelse(test = img_orientation == "portrait", 
                               yes = 0.3, no = 0.05),
                xright = ifelse(test = img_orientation == "portrait", 
                                yes = 0.7, no = 0.95),
                ybottom = ifelse(test = img_orientation == "portrait", 
                                 yes = 0.3, no = 0.4),
                ytop = ifelse(test = img_orientation == "portrait", 
                              yes = 0.7, no = 0.6))
  }
 
  
  box()
  
  for (igear in sort(available_gears, decreasing = T) ) {
    
    par(mar = c(0, 0, 0, 0))
    
    ## Load full VAST result
    fit <- readRDS(paste0("results/chukchi_", igear, "/vast_fits/", 
                          ispp_name, "/fit_full.rds"))
    
    
    ## Years that gear was used
    years <- list("otter" = data.frame(idx = c(1, 23),
                                       year = c(1990, 2012)),
                  "beam" = data.frame(idx = c(1, 6, 8),
                                      year = c(2012, 2017, 2019)))[[igear]]
    
    temp_density <- get(paste0("ms_dens_", igear))[, ispp_name, ]
    density_cuts <- round(c(0, 1, 
                            quantile(x = temp_density[temp_density > 1],
                                     probs = c(0.25, 0.5, 0.75)),
                            max(temp_density) + 1), 1)
    assign(x = paste0("density_cuts_", igear), value = density_cuts)
    
    for (iyear in 1:nrow(years)){
      plot_this <- data.frame(X1 = temp_density[, iyear])
      
      plot_this$X1 <- as.numeric(base::cut(x = plot_this$X1, 
                                           breaks = density_cuts, 
                                           labels = 1:(length(density_cuts) - 1),
                                           include.lowest = TRUE,
                                           right = FALSE))
      chukchi_shp <- sp::SpatialPointsDataFrame(coords = grid_pts@coords, 
                                                proj4string = aea_crs,
                                                data = plot_this)
      chukchi_ras <- raster::raster(x = chukchi_shp, resolution = 3700)
      chukchi_ras <- raster::rasterize(x = chukchi_shp,
                                       y = chukchi_ras, 
                                       field = "X1")
      
      image(chukchi_ras, asp = 1, col = get(paste0("col_", igear)), 
            axes = F, ann = F)
      box()
      
      ## Plot Observed densities
      assign(x = paste0("obs_dens_", igear), 
             value = subset(x = fit$data_frame[fit$data_list$PredTF_i == 0, ], 
                            subset = t_i == years$year[iyear]))
      latlon_locs <- 
        sp::SpatialPoints(coords = get(paste0("obs_dens_", 
                                              igear))[, c("Lon_i", "Lat_i")],
                          proj4string = latlon_crs)
      aea_locs <- sp::spTransform(x = latlon_locs, CRSobj = aea_crs)
      obs_dens <- get(paste0("obs_dens_", igear))$b_i
      points(aea_locs[obs_dens == 0], pch = 3, cex = 0.5)
      points(aea_locs[obs_dens < 1 & obs_dens > 0], cex = 0.1)
      points(aea_locs[obs_dens > 0], cex = log10(obs_dens[obs_dens > 0]) / 2)
      
      ## Plot land
      plot(ak_land_cropped, add = T, border = F, col = "tan")
      
      ## Plot species/gear label
      text(x = extent(chukchi_ras)[1] + diff(extent(chukchi_ras)[1:2])*0.75,
           y = extent(chukchi_ras)[3] + diff(extent(chukchi_ras)[3:4])*0.3,
           labels = paste0(igear, " trawl\n", years$year[iyear]),
           cex = 1)
    } 
    
    ## Index
    index <- read.csv(paste0("results/chukchi_", igear, "/vast_fits/", 
                             ispp_name, "/Table_for_SS3.csv" ))[years$idx, ]
    
    index_est <- rbind(index$Estimate_metric_tons - index$SD_mt,
                       index$Estimate_metric_tons,
                       index$Estimate_metric_tons + index$SD_mt) 
    
    scale_est <- ifelse(test = max(log10(abs(index_est))) > 6, 
                        yes = "million",
                        no = "thousand")
    
    index_est <- index_est / c("thousand" = 1e3, "million" = 1e6)[scale_est]
    
    year_bounds <- list("otter" = c(1985, 2017),
                        "beam" = c(2010, 2021))[[igear]]
    
    par(mar = c(4, 4, 2, 0.5))
    plot(x = index$Year,
         y = index_est[2, ],
         axes = F,
         yaxs = "i", yaxt = "n",
         xlim = year_bounds, 
         ylim = c(0, 1.15 * max(index_est[3, ])), las = 1, 
         pch = 16, col = get(paste0("col_", igear))[8],
         xlab = "Year", ylab = "")
    axis(side = 2, las = 1, cex.axis = 0.75)
    mtext(side = 2, text = paste0("Index (", scale_est, "\nmetric tons)"), 
          line = 2, cex = 0.5)
    segments(x0 = index$Year, x1 = index$Year,
             y0 = index_est[1, ], y1 = index_est[3, ],
             col = get(paste0("col_", igear))[7])
    text(x = index$Year, y = index_est[3, ], xpd = NA,
         labels = paste("CV =", round(index$SD_log, 2)),
         pos = 3, cex = 0.7, col = get(paste0("col_", igear))[8])
    # axis(side = 2, las = 1)
    axis(side = 1, at = c(0, years$year), labels = c(0, years$year), las = 2)
    box(which = "figure")
    
    ###################################
    ## Calculate DHarma Residuals
    ###################################
    dyn.load(paste0("results/chukchi_", igear, "/vast_fits/", 
                    ispp_name, "/VAST_v12_0_0.dll"))
    
    dharmaRes = summary( fit, what = "residuals", working_dir = NA )
    dyn.unload(paste0("results/chukchi_", igear, "/vast_fits/", 
                      ispp_name, "/VAST_v12_0_0.dll"))
    
    par(mar = c(3.5, 4.5, 1.25, 1))
    
    ###################################
    ## QQ Plot
    ###################################
    gap::qqunif(dharmaRes$scaledResiduals, 
                pch = 2, 
                bty = "n",
                logscale = F, 
                col = "black", 
                cex = 0.2,
                cex.main = 1, 
                las = 1,
                ann = F, 
                cex.axis = 1)
    
    mtext(side = 1, line = 2, text = "Expected", cex = 0.7)
    mtext(side = 2, line = 2.5, text = "Observed", cex = 0.7)
    
    box()
    box(which = "figure")
    
    ## Legend
    par(mar = c(0,0,0,0))
    plot(1, type = "n", xlim = c(0, 100), ylim = c(0, 100), ann = F, axes = F)

    if (igear == "otter")
      plotrix::color.legend(
        xl = 10, xr = 20,
        yb = 5, yt = 85,
        legend = c("<1", round(density_cuts_otter)[-(1:2)]),
        rect.col = get(paste0("col_otter")),
        gradient = "y",
        align = "rb",
        cex = 0.4)

    if (igear == "beam")
      plotrix::color.legend(
        xl = 10, xr = 20,
        yb = 5, yt = 85,
        legend = c("<1", round(density_cuts_beam)[-(1:2)]),
        rect.col = get(paste0("col_beam")),
        gradient = "y",
        align = "rb",
        cex = 0.4)

    points(x = rep(70, 6),
           y = seq(from = 10, to = 80, length = 6),
           pch = c(3, rep(1, 5)),
           cex = c(0.5, 0.1, 1:4 / 2))
    text(x = rep(70, 6),
         y = seq(from = 10, to = 80, length = 6) - 5,
         labels = c(0, "< 1", 10, 100, 1000, 10000),
         cex = 0.65)

    mtext(side = 3, line = -1.75, cex = 0.55,
          text = paste0(igear, " trawl\n(kg/km2)"))
    box()
  }
  

  dev.off()
  
}
