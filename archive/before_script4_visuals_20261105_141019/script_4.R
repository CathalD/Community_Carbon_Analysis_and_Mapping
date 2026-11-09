
# ==========================

# Part 3: prepare the prior and compare mapped stocks

# ===============
library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/grid_helper.R")

# 1. Inputs and the shared analysis grid
output_folder <- "outputs/script4"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask

# Transfer external values to the shared cell centres once.
align_to_grid <- function(x) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    x <- project(x, template, method = "near")
  }
  mask(x, analysis_mask)
}

# 2. Prepare paired prior stock and uncertainty
prepared_priors <- list()
for (depth in intersect(names(prior_mean_files), names(prior_width_files))) {
  source_mean <- rast(prior_mean_files[[depth]])
  source_width <- rast(prior_width_files[[depth]])
  compareGeom(source_mean, source_width)
  pair <- align_to_grid(c(source_mean, source_width))
  names(pair) <- c("prior_mean", "interval_width_90")
  prior_sd <- pair$interval_width_90 / (2 * qnorm(0.95))
  usable <- is.finite(pair$prior_mean) & pair$prior_mean > 0 &
    is.finite(prior_sd) & prior_sd > 0
  pair <- mask(pair, usable, maskvalues = 0)
  prior_sd <- mask(prior_sd, usable, maskvalues = 0)
  log_variance <- log1p((prior_sd / pair$prior_mean)^2)
  log_mean <- log(pair$prior_mean) - log_variance / 2
  prior <- c(pair$prior_mean, prior_sd, log_mean, log_variance)
  names(prior) <- c("prior_mean", "prior_sd", "log_mean", "log_variance")
  writeRaster(prior, file.path(output_folder, paste0("prior_moments_0_", depth, "cm.tif")),
              datatype = "FLT8S", overwrite = TRUE)
  # Keep the familiar filenames; the width layer retains its original units.
  writeRaster(pair$prior_mean, file.path(output_folder, basename(prior_mean_files[[depth]])),
              overwrite = TRUE)
  writeRaster(pair$interval_width_90,
              file.path(output_folder, basename(prior_width_files[[depth]])), overwrite = TRUE)
  prepared_priors[[depth]] <- prior
  cat("\nPrior 0–", depth, "cm: Q95-Q05 -> approximate SD -> lognormal moments.\n", sep = "")
  print(global(prior, c("min", "max", "notNA"), na.rm = TRUE))
  png(file.path(output_folder, paste0("prior_0_", depth, "cm.png")),
      width = 1400, height = 700, res = 150)
  plot(prior[[1:2]], main = c("Prior mean (kg C/m²)", "Approximate prior SD (kg C/m²)"))
  dev.off()
  plot(prior[[1:2]], main = c("Prior mean", "Approximate prior SD"))
}

# 3. Compare prediction maps at matching depths and on shared valid cells
# Script 3's native analysis is unchanged. Its selected raster is transferred
# to the common grid here, rather than changing its covariates or fitted models.
for (depth in mapping_depths_cm) {
  suffix <- paste0("0_", depth, "cm")
  files <- c(
    Area_average = file.path("outputs/script2", paste0("area_based_", suffix, ".tif")),
    IDW = file.path("outputs/script2", paste0("idw_", suffix, ".tif"))
  )
  ok_file <- file.path("outputs/script2", paste0("ok_", suffix, ".tif"))
  if (file.exists(ok_file)) files <- c(files, Ordinary_kriging = ok_file)
  maps <- lapply(files, function(path) {
    x <- rast(path)
    compareGeom(x, template)
    mask(x, analysis_mask)
  })
  regression_folder <- file.path("outputs/script3", paste0("stock_0_", depth))
  regression_file <- file.path(regression_folder, "selected_carbon_stock.tif")
  info_file <- file.path(regression_folder, "selected_model.csv")
  if (file.exists(regression_file) && file.exists(info_file)) {
    info <- read_csv(info_file, show_col_types = FALSE)
    regression <- align_to_grid(rast(regression_file))
    maps[[paste0("Regression_", info$model[1])]] <- regression
    writeRaster(regression,
                file.path(output_folder, paste0("regression_selected_", suffix, ".tif")),
                overwrite = TRUE)
  }
  if (as.character(depth) %in% names(prepared_priors)) {
    maps$Prior <- prepared_priors[[as.character(depth)]][["prior_mean"]]
  }
  comparison <- do.call(c, unname(maps))
  names(comparison) <- names(maps)
  shared_cells <- app(is.finite(comparison), "sum") == nlyr(comparison)
  comparison <- mask(comparison, shared_cells, maskvalues = 0)
  map_summary <- global(comparison, c("mean", "min", "max", "notNA"), na.rm = TRUE)
  cat("\nMap comparison:", suffix, "(map differences describe agreement, not accuracy)\n")
  print(map_summary)
  writeRaster(comparison, file.path(output_folder, paste0("comparison_0_", depth, ".tif")),
              overwrite = TRUE)
  limits <- range(minmax(comparison), finite = TRUE)
  png(file.path(output_folder, paste0("maps_0_", depth, ".png")),
      width = 500 * nlyr(comparison), height = 650, res = 150)
  plot(comparison, range = limits)
  dev.off()
  plot(comparison, range = limits)
}
cat("\nPrepared prior and comparison rasters:", output_folder, "\n")
