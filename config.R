# Shared settings for the workshop. Run scripts from the project folder.

analysis_crs <- "EPSG:3978"
coordinate_crs <- "EPSG:4326"
cell_size_m <- 1000
grid_folder <- "data/grid"

workbook_path <- "data/field/field_cores.xlsm"
mapping_file <- "outputs/script1/mapping_data.csv"
aoi_file <- "data/boundary/aoi.geojson"
elevation_file <- "data/covariates/topography.tif"
water_file <- "data/covariates/water.tif"

# Use Script 1 stocks unchanged, capped at the sampled maximum.
# Field profiles are assumed complete at that maximum; do not exclude shallow cores.
mapping_depths_cm <- c(30, 50, 100)

# Script 3 retains its existing response scale, models and native raster analysis.
regression_stock_variable <- "stock_0_30"
water_threshold <- 50

# A matching stock AND uncertainty prior must exist for the selected depth.
bayesian_depth_cm <- 100
prior_mean_files <- c("100" = "data/priors/mcmaster_0_100cm.tif")
prior_width_files <- c("100" = "data/priors/mcmaster_uncertainty_0_100cm.tif")

grid_file <- file.path(grid_folder, "template_1km.tif")
mask_file <- file.path(grid_folder, "aoi_mask_1km.tif")
cell_id_file <- file.path(grid_folder, "cell_ids_1km.tif")
projected_aoi_file <- file.path(grid_folder, "aoi_projected.gpkg")
