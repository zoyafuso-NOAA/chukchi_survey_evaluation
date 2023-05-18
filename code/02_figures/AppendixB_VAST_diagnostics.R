##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Appendix B: VAST output
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot predicted density, index, qq plots for each species/gear
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import relevant libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(terra)
library(RColorBrewer)
library(gap)
library(plotrix)
library(FishStatsUtils)
library(jpeg)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data and constants ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(file = "data/survey_opt_data/optimization_data.RData")
AK_land <- terra::vect(x = "data/spatial_data/land_shapefiles/AKland.shp")
grid_pts <- terra::vect(x = "data/survey_opt_data/grid_pts.shp")
cropped_extent <- terra::ext(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- terra::crop(x = AK_land, y = cropped_extent)
good_species <- read.csv(file = "results/good_species.csv")

## Plotting colors
col_otter <- RColorBrewer::brewer.pal(n = 9, name = "Blues")
col_beam <- RColorBrewer::brewer.pal(n = 9, name = "Greens")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Result directory ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
if(!dir.exists(paths = "figures/Appendix_B/"))
  dir.create(path = "figures/Appendix_B/")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Collate density ----
##    predictions for each gear and species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (igear in c("otter", "beam")) {
  
  ## Years that gear was used
  years <- list("otter" = data.frame(idx = c(1, 23),
                                     year = c(1990, 2012)),
                "beam" = data.frame(idx = c(1, 6, 8),
                                    year = c(2012, 2017, 2019)))[[igear]]
  
  ## Species that gear sampled
  species <- read.csv(file = "results/good_species.csv")
  species_idx <- species[, igear]
  species <- species[species_idx, ]$taxon
  
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (irow in 1:nrow(x = good_species)) {

  ## Taxon constants
  ispp_name <- good_species$taxon[irow]
  ispp_name_plot <- good_species$taxon_plot[irow]
  available_gears <- unlist(subset(x = good_species,
                                   subset = taxon == ispp_name,
                                   select = c(otter, beam)))
  available_gears <- sort(names(available_gears)[available_gears])
  n_gears <- length(available_gears)
  
  ## Open plot device
  png(filename = paste0("figures/Appendix_B/VAST_output_", ispp_name, ".png"),
      width = 190, units = "mm", res = 500, family = "serif",
      height = switch(paste0(available_gears, collapse = ""),
                      "beam" =  40, "otter" = 45,
                      "beamotter" = 75))
  
  ## Setup plot layout 
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
                          "beam" = c(0.75, 1, 1, 1, 0.5, 1, 1),
                          "otter" = c(1, 1, 1, 0.5, 1, 1))
  
  layout(mat = plot_layout_mat,
         widths = layout_widths)
  par(mar = c(0, 0, 0, 0))
  
  ## Taxon name 
  plot(1, axes = F, ann = F, pch = 16, type = "n",
       xlim = c(0, 1), ylim = c(0, 1))
  mtext(side = 3,
        text = gsub(x = ispp_name_plot,
                    pattern = "\\n",
                    replacement = "\n", fixed = TRUE),
        line = -3)
  
  ## Taxon image
  if(file.exists(paste0("data/taxon_images/", ispp_name, ".jpg"))) {
    
    xleft <- good_species[good_species$taxon == ispp_name, "xleft"]
    xright <- good_species[good_species$taxon == ispp_name, "xright"]
    ybottom <- good_species[good_species$taxon == ispp_name, "ybottom"]
    ytop <- good_species[good_species$taxon == ispp_name, "ytop"]
    
    fish_img <- jpeg::readJPEG(source = paste0("data/taxon_images/", 
                                               ispp_name, ".jpg"))
    img_orientation <- ifelse(test = dim(fish_img)[1] > dim(fish_img)[2],
                              yes = "portrait", no = "landscape")
    rasterImage(image = fish_img, 
                xleft = xleft,
                xright = xright,
                ybottom = ybottom,
                ytop = ytop)
  }
  
  box()
  
  for (igear in sort(available_gears, decreasing = T) ) {
    
    par(mar = c(0, 0, 0, 0))
    
    ## Load full VAST result
    fit <- readRDS(file = paste0("results/chukchi_", igear, "/vast_fits/", 
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
      
      plot_this$X1 <- 
        as.numeric(base::cut(x = plot_this$X1, 
                             breaks = density_cuts, 
                             labels = 1:(length(density_cuts) - 1),
                             include.lowest = TRUE,
                             right = FALSE))
      
      chukchi_shp <- terra::vect(x = cbind(chukchi_sea_grid,
                                           plot_this), 
                                 geom = c("Lon", "Lat"), 
                                 crs = latlon_crs)
      chukchi_shp <- terra::project(x = chukchi_shp, aea_crs)
      
      chukchi_ras <- terra::rast(x = chukchi_shp, resolution = 4000)
      chukchi_ras <- terra::rasterize(x = chukchi_shp,
                                      y = chukchi_ras, 
                                      field = "X1")

      image(x = chukchi_ras, asp = 1, col = get(paste0("col_", igear)), 
            axes = F, ann = F)
      box()
      
      ## Plot Observed densities
      assign(x = paste0("obs_dens_", igear), 
             value = subset(x = fit$data_frame[fit$data_list$PredTF_i == 0, ], 
                            subset = t_i == years$year[iyear]))
      latlon_locs <- terra::vect(x = get(paste0("obs_dens_", igear)),
                                 geom = c("Lon_i", "Lat_i"),
                                 crs = latlon_crs)
      aea_locs <- terra::project(x = latlon_locs, aea_crs)

      obs_dens <- get(paste0("obs_dens_", igear))$b_i
      points(aea_locs[obs_dens == 0], pch = 3, cex = 0.5)
      points(aea_locs[obs_dens < 1 & obs_dens > 0], cex = 0.1, pch = 1)
      points(aea_locs[obs_dens > 0], 
             cex = log10(obs_dens[obs_dens > 0]) / 2, pch = 1)
      
      ## Plot land
      plot(ak_land_cropped, add = T, border = F, col = "tan")
      
      ## Plot species/gear label
      text(x = cropped_extent[1] + diff(cropped_extent[1:2])*0.75,
           y = cropped_extent[3] + diff(cropped_extent[3:4])*0.3,
           labels = paste0(igear, " trawl\n", years$year[iyear]),
           cex = 1)
    } 
    
    ## Index
    index <- read.csv(file = paste0("results/chukchi_", igear, "/vast_fits/", 
                                    ispp_name, 
                                    "/diagnostics/Index.csv" ))[years$idx, ]
    
    index_est <- rbind(index$Estimate - index$Std..Error.for.Estimate,
                       index$Estimate,
                       index$Estimate + index$Std..Error.for.Estimate) * 1e-3
    
    scale_est <- ifelse(test = max(log10(abs(index_est / 1000))) > 6, 
                        yes = "million",
                        no = "thousand")
    
    index_est <- index_est / c("thousand" = 1e3, "million" = 1e6)[scale_est]
    
    year_bounds <- list("otter" = c(1980, 2017),
                        "beam" = c(2010, 2021))[[igear]]
    
    par(mar = c(4, 3.5, 2, 0.5))
    plot(x = index$Time,
         y = index_est[2, ],
         axes = F,
         yaxs = "i", yaxt = "n",
         xlim = year_bounds, 
         ylim = c(0, 1.15 * max(index_est[3, ])), las = 1, 
         pch = 16, col = get(paste0("col_", igear))[8],
         xlab = "Year", ylab = "")
    axis(side = 2, las = 1, cex.axis = 0.75)
    mtext(side = 2, text = paste0("Index (", scale_est, " metric tons)"), 
          line = 2.5, cex = 0.5)
    segments(x0 = index$Time, x1 = index$Time,
             y0 = index_est[1, ], y1 = index_est[3, ],
             col = get(paste0("col_", igear))[7])
    text(x = index$Time, y = index_est[3, ], xpd = NA,
         labels = paste("CV =", round(index$Std..Error.for.ln.Estimate., 2)),
         pos = 3, cex = 0.7, col = get(paste0("col_", igear))[8])
    axis(side = 1, at = c(0, years$year), labels = c(0, years$year), las = 2)
    box(which = "figure")
    
    ## Calculate DHarma Residuals
    load(file = paste0("results/chukchi_", igear, "/vast_fits/",
                       ispp_name, "/diagnostics/diagnostics.RData"))
    dharmaRes <- diagnostics$dharmaRes
    
    ## QQ Plot
    par(mar = c(3.5, 3.5, 1, 1))
    gap::qqunif(dharmaRes$scaledResiduals, 
                pch = 2, 
                bty = "n",
                logscale = F, 
                col = "black", 
                cex = 0.2,
                cex.main = 1, 
                las = 2,
                ann = F, 
                cex.axis = 0.7)
    
    mtext(side = 1, line = 2, text = "Expected", cex = 0.7)
    mtext(side = 2, line = 2, text = "Observed", cex = 0.7)
    
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
