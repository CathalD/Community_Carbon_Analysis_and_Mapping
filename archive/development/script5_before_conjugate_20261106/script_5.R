# ===================================

# Part 2.3: Bayesian carbon updates

# ========================

library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/grid_helper.R")


# Set-up out settings here
bayesian_depth_cm = 100 # We will start with 100cm, the we can iterate
target_depth_cm <- bayesian_depth_cm
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
range_m <- 2000 # Can be adjusted - Note to self, should run a maximum likelihood


observation_log_sd <- 0 # Set observation error to zero for now
chunk_size <- 5000L 
map_methods <- c("area_based", "idw") # If any outputs work as long as there is a mean and sd
output_folder <- file.path("outputs/script5", paste0("stock_0_", target_depth_cm))
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# 2. Read the shared grid and prior already prepared by Script 4
shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask
prior <- rast(file.path("outputs/script4", paste0("prior_moments_0_", target_depth_cm, "cm.tif")))
compareGeom(prior, template)
prior <- mask(prior, analysis_mask)
prior_values <- values(prior, mat = TRUE)
eligible_cells <- which(apply(is.finite(prior_values), 1, all) &
                          prior_values[, "prior_mean"] > 0 &
                          prior_values[, "log_variance"] > 0)
grid_xy <- xyFromCell(template, eligible_cells)
grid_mean <- prior_values[eligible_cells, "prior_mean"]
grid_sd <- prior_values[eligible_cells, "prior_sd"]
grid_log_mean <- prior_values[eligible_cells, "log_mean"]
grid_log_var <- prior_values[eligible_cells, "log_variance"]

raw_cores <- read_csv(mapping_file, show_col_types = FALSE)
cores <- data.frame(
  core_id = raw_cores[["Core ID"]], plot_id = raw_cores[["Plot ID"]],
  longitude = raw_cores$Longitude, latitude = raw_cores$Latitude,
  stock = raw_cores[[stock_column]]
)
cores <- cores[is.finite(cores$stock) & cores$stock > 0, ]
core_points <- project(vect(cores, geom = c("longitude", "latitude"), crs = coordinate_crs),
                       crs(template))
core_xy <- crds(core_points)

# Lookup the source prior at the actual core coordinates, not at moved centres.
# The final predictions still use the shared analysis grid.
depth_key <- as.character(target_depth_cm)
native_prior <- c(rast(prior_mean_files[[depth_key]]),
                  rast(prior_width_files[[depth_key]]) / (2 * qnorm(0.95)))
names(native_prior) <- c("prior_mean", "prior_sd")
at_cores <- extract(native_prior, project(core_points, crs(native_prior)), ID = FALSE)
core_log_var <- log1p((at_cores$prior_sd / at_cores$prior_mean)^2)
core_log_mean <- log(at_cores$prior_mean) - core_log_var / 2
cores$prior_mean <- at_cores$prior_mean
cores$prior_sd <- at_cores$prior_sd
cat("\nTarget:", stock_column, "Cores:", nrow(cores),
    "Plots:", length(unique(cores$plot_id)), "\n")
cat("Raw-core range:", range_m, "m; observation log SD:", observation_log_sd, "\n")
print(cores)


