##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Figure X: Plot STRS optimized solutions
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot STRS optimization solutions with different allocations
##                     based on species tradeoffs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import relevant libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(rgdal)
library(sp)
library(raster)
library(RColorBrewer)
library(graticule)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##   Crop the Alaska shapefile to only include NW AK
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spp_settings <- read.csv(file = "results/good_species.csv", 
                         stringsAsFactors = F)
load(file = "data/survey_opt_data/optimization_data.RData")
ak_land <- rgdal::readOGR(dsn = "data/spatial_data/land_shapefiles/AKland.shp")
cropped_extent <- extent(grid_pts)
cropped_extent[2] <- cropped_extent[2] + 10000
ak_land_cropped <- raster::crop(x = ak_land, y = cropped_extent)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import STRS solutions optimized for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_strata <- c("beam" = 4, "otter" = 3)
load(file = paste0("results/chukchi_otter/survey_opt/Str_", 
                   n_strata["otter"], "/result_list.RData"))
assign(x = "otter_result_list", value = result_list)
load(file = paste0("results/chukchi_otter/survey_opt/Str_", 
                   n_strata["otter"], "/allocations.RData"))
assign(x = "otter_allocations", value = ss_sample_allocations)

load(file = paste0("results/chukchi_beam/survey_opt/Str_", 
                   n_strata["beam"], "/result_list.RData"))
assign(x = "beam_result_list", value = result_list)
load(file = paste0("results/chukchi_beam/survey_opt/Str_", 
                   n_strata["beam"], "/allocations.RData"))
assign(x = "beam_allocations", value = ss_sample_allocations)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set up graticules ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
aea_grat <- graticule::graticule(
  lons = seq(from = -180, to = -150, by = 5),
  lats = seq(from = 65, to = 75, by = 2.5), 
  proj = aea_crs)
aea_grat_labs <- graticule::graticule_labels(
  lons = seq(from = -180, to = -150, by = 10),
  lats = seq(from = 70, to = 72.5, by = 2.5),
  xline = -171, yline = 72.5,
  proj = aea_crs)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Which species were represented by both gears or just one gear?
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# two_spp_idx <- which(x = rowSums(spp_settings[, c("beam", "otter")]) == 2)

# two_spp <- spp_settings$taxon[two_spp_idx]
# n_two_spp <- length(two_spp)
# two_spp_title <- spp_settings$taxon_plot[two_spp_idx]
# two_spp_title <- gsub(x = two_spp_title, pattern = "\\n", replacement = "\n",
#                       fixed = TRUE)

otter_beam_idx <- which(spp_settings$otter == T & spp_settings$beam == T)
otter_spp_idx <- which(spp_settings$otter == T & spp_settings$beam == F)
beam_spp_idx <- which(spp_settings$otter == F & spp_settings$beam == T)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot ----
##   There are 12 species represented by both gears. 
##   There are 5 species just represented by the otter trawl: Alaska place, 
##         eelpouts, P. herring, walleye pollock, yellowfin sole.
##   There are 3 species just represented by the beam trawl: pricklebacks,
##         bryozoans, and bivalves.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
{
  ## Open png device
  png(filename = paste0("figures/FigureX_ss_strs_sols.png"),
      units = "mm", width = 190, height = 155, res = 500)
  
  layout(mat = matrix(data = c(1, 2, 31, 3, 4, 31, 5, 6, 31, 23,27,
                               7, 8, 31, 9,10, 31, 11,12,31, 24,28,
                               13,14,31, 15,16,31, 17,18,31, 25,29,
                               31,31,31, 19,20,31, 21,22,31,  30,26),
                      nrow = 4, byrow = TRUE),
         widths = c(1,1,0.2, 1,1,0.2, 1,1,0.2, 1,1))
  
  par(mar = c(0, 0, 2.5, 0),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = "serif")
  
  for (ispp in c(otter_beam_idx, 
                 otter_spp_idx, 
                 beam_spp_idx)) { ## Loop over species -- start
    for (igear in c("otter", "beam")) { ## Loop over gears -- start
      
      if (igear == "otter" & ispp %in% beam_spp_idx) next
      if (igear == "beam" & ispp %in% otter_spp_idx) next
      
      ## temporarily name result list and allocations
      result_list <- get(paste0(igear, "_result_list"))
      ss_sample_allocations <- get(paste0(igear, "_allocations"))
      
      ## Create solution shape object
      solution <- sp::SpatialPointsDataFrame(
        coords = grid_pts@coords,
        data = data.frame(Str_no = result_list$sol_by_cell), 
        proj4string = aea_crs )
      solution_ras <- raster::raster(x = solution,
                                     resolution = 4000)
      solution_ras <- raster::rasterize(x = solution,
                                        y = solution_ras,
                                        field = "Str_no")
      strata_cols <- RColorBrewer::brewer.pal(n = n_strata[igear], 
                                              name = "Set2") 
      
      ## Subset 100-sample allocations for plotting purposes
      ss_solutions <- subset(x = ss_sample_allocations,
                             subset = n == 100 & 
                               species == spp_settings$taxon[ispp])
      
      ## Take a random sample based on the STRS design optimized for ispp
      sample_vec <- c()
      for (istratum in 1:n_strata[igear]) { ## Loop over strata -- start
        sample_vec <- 
          c(sample_vec,
            sample(x = which(x = result_list$sol_by_cell == istratum),
                   size = ss_solutions[, paste0("Str_", istratum)]))
      } ## Loop over strata -- end
      
      ## Turn randomly drawn stations into a spatial object
      samp_pts <- sp::SpatialPoints(coords = grid_pts@coords[sample_vec, ], 
                                    proj4string = aea_crs)
      
      ## Plot solution
      image(x = solution_ras, 
            col = strata_cols, asp = 1, axes = F, ann = F)
      
      ## Plot sampled stations
      points(samp_pts, pch = 16, cex = 0.75)
      
      ## Plot graticules
      plot(aea_grat, add = T, lwd = 0.5)
      
      ## Plot land
      plot(ak_land_cropped, col = "tan", border = F, add = TRUE)
      
      ## Plot species title
      if (ispp %in% otter_beam_idx & igear == "otter")
        text(x = extent(solution_ras)[1] + 
               diff(extent(solution_ras)[1:2]) * 1,
             y = extent(solution_ras)[3] + 
               diff(extent(solution_ras)[3:4]) * 1.225,
             labels = gsub(pattern = '\\n', 
                           x = spp_settings$taxon_plot[ispp], 
                           replacement = "\n",
                           fixed = TRUE), 
             pos = 1, xpd = NA, cex = 1.25, font = 2)
      
      if (ispp %in% c(otter_spp_idx, beam_spp_idx))
        mtext(side = 3, line = 0.75,
              text = spp_settings$taxon_plot[ispp],
              cex = 0.75, font = 2)
      
      ## Plot graticule labels
      text(subset(x = aea_grat_labs, 
                  subset = lab %in% paste0(c(170,160), "*degree*W")),
           lab = parse(text = aea_grat_labs$lab[aea_grat_labs$lab %in% 
                                                  paste0(c(170, 160), 
                                                         "*degree*W")]),
           pos = 3, cex = 0.45, xpd = NA)
      text(subset(aea_grat_labs, !aea_grat_labs$islon), 
           lab = parse(text = aea_grat_labs$lab[!aea_grat_labs$islon]),
           pos = 1, cex = 0.45, xpd = NA)
      
      ## Scale bar
      segments(x0 = extent(x = solution_ras)[1] + 
                 diff(x = extent(x = solution_ras)[1:2]) * 0.75, 
               x1 = extent(x = solution_ras)[1] + 
                 diff(x = extent(x = solution_ras)[1:2]) * 0.75 + 100000,
               y0 = extent(x = solution_ras)[3] + 
                 diff(x = extent(x = solution_ras)[3:4]) * 0.1,
               lwd = 2)
      text(x = extent(x = solution_ras)[1] + 
             diff(x = extent(x = solution_ras)[1:2]) * 0.83,
           y = extent(x = solution_ras)[3] + 
             diff(x = extent(x = solution_ras)[3:4]) * 0.11,
           labels = "100 km", 
           pos = 1, cex = 0.6)
      
      ## Sampling Density Legend
      legend(x = extent(x = solution_ras)[1] + 
               diff(x = extent(x = solution_ras)[1:2]) * 0.45,
             y = extent(x = solution_ras)[3] + 
               diff(x = extent(x = solution_ras)[3:4]) * 0.525,
             legend = rev(x = paste0(
               "Str ", 1:n_strata[igear], ": ", 
               ss_solutions[, paste0("Str_", 1:n_strata[igear])]/100)),
             title = paste0(igear, " trawl"),
             fill = rev(x = strata_cols), cex = 0.6, bty = "n")
      
    } ## Loop over gears -- end
  } ## Loop over species -- end
  dev.off()
}
