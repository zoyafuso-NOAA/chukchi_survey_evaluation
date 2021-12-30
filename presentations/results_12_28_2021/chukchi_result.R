source("modified_functions/plot_survey_opt_map.R")

istrata = 3
chukchi_grid <- read.csv("data/spatial_data/BS_Chukchi_extrapolation_grids/ChukchiThorsonGrid.csv")
load(paste0("results/chukchi_beam/survey_opt/MS/Str_", 
            istrata, "/result_list.RData"))
load(paste0("results/chukchi_beam/survey_opt/MS/Str_", 
            istrata, "/allocations.RData"))
load("data/survey_opt_data/optimization_data.RData")

plot_survey_opt_map(file_name = paste0("presentations/results_12_28_2021/",
                                       "base_solution.png"),
                    grid_object =  grid_pts,
                    sol_by_cell = result_list$sol_by_cell, 
                    allocations = as.numeric(ms_sample_allocations[ms_sample_allocations$n == 100, paste0("Str_ ", 1:istrata)]),
                    draw_stations = FALSE)

## loop over species
for (ispp in 1:n_spp_beam) {
  temp_sol <- as.integer(subset(ss_sample_allocations, species == spp_list_beam[ispp] & n == 70, select = paste("Str_", 1:istrata)))
  
  plot_survey_opt_map(file_name = paste0("presentations/results_12_28_2021/",
                                         spp_list_beam[ispp],
                                         "_ss_solution.png"),
                      grid_object =  grid_pts,
                      sol_by_cell = result_list$sol_by_cell, 
                      allocations = temp_sol,
                      draw_stations = TRUE)
}

temp_sol <- as.integer(ms_sample_allocations[ms_sample_allocations$n == 70, paste0("Str_ ", 1:istrata)])

plot_survey_opt_map(file_name = paste0("presentations/results_12_28_2021/",
                                       "ms_solution.png"),
                    grid_object =  grid_pts,
                    sol_by_cell = result_list$sol_by_cell, 
                    allocations = temp_sol,
                    draw_stations = TRUE)
tapply(X = chukchi_grid$Shape_Area * 1e-6, INDEX = result_list$sol_by_cell, FUN = sum) / sum(chukchi_grid$Shape_Area * 1e-6)
temp_sol / sum(temp_sol)
