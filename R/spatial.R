# Small, illustrative spatial examples. No raster downloads or external services.
# All stock comparisons use soil organic carbon at 0–30 cm, in kg C/m2.

.spatial_packages <- function() {
  for (p in c("sf", "terra", "ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, " before this step.")
  }
}

.spatial_metadata <- function(out, method, formula, packages, assumptions) {
  utils::write.csv(data.frame(method = method, formula = formula,
    packages = packages, assumptions = assumptions),
    file.path(out, "method_metadata.csv"), row.names = FALSE)
}

.read_spatial_inputs <- function(config, core_path = "outputs/part1/core_stocks.csv") {
  .spatial_packages()
  if (!isTRUE(config$illustrative)) {
    stop("These maps require illustrative = TRUE; their ecosystems and uncertainty are assumptions.")
  }
  if (!file.exists(core_path)) stop("Run scripts/01_field.R first.")
  if (is.null(config$aoi_path) || !file.exists(config$aoi_path)) stop("AOI GeoJSON is missing.")
  cores <- utils::read.csv(core_path, stringsAsFactors = FALSE)
  needed <- c("core_id", "ecosystem", "group", "latitude", "longitude", "stock_0_30_kg_m2")
  if (!all(needed %in% names(cores))) stop("Core stock table lacks spatial input columns.")
  if (anyDuplicated(cores$core_id)) stop("Core IDs must be unique.")
  if (any(!is.finite(cores$latitude) | !is.finite(cores$longitude) |
    !is.finite(cores$stock_0_30_kg_m2))) stop("Each mapped core needs coordinates and a complete 0–30 cm stock.")
  if (any(abs(cores$latitude) > 90 | abs(cores$longitude) > 180)) stop("Invalid latitude or longitude.")
  aoi <- sf::st_read(config$aoi_path, quiet = TRUE)
  if (is.na(sf::st_crs(aoi))) stop("The AOI needs a coordinate reference system.")
  aoi <- sf::st_transform(sf::st_make_valid(aoi), 4326)
  bb <- sf::st_bbox(aoi)
  local_crs <- paste0("+proj=laea +lat_0=", (bb["ymin"] + bb["ymax"]) / 2,
    " +lon_0=", (bb["xmin"] + bb["xmax"]) / 2, " +datum=WGS84 +units=m +no_defs")
  aoi <- sf::st_sf(area_id = "Illustrative study area", geometry = sf::st_union(sf::st_transform(aoi, local_crs)))
  points <- sf::st_as_sf(cores, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  points <- sf::st_transform(points, local_crs)
  if (!all(lengths(sf::st_intersects(points, aoi)) > 0)) stop("A mapped core is outside the AOI.")
  list(aoi = aoi, points = points, cores = cores)
}

.area_grid <- function(aoi, cellsize_m = 500) {
  if (length(cellsize_m) != 1 || !is.finite(cellsize_m) || cellsize_m <= 0) stop("Invalid grid resolution.")
  grid <- sf::st_make_grid(aoi, cellsize = cellsize_m, square = TRUE)
  full <- sf::st_sf(cell_id = seq_along(grid), geometry = grid)
  centres <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(full)))
  full$x_m <- centres[, 1]; full$y_m <- centres[, 2]
  full <- full[lengths(sf::st_intersects(full, aoi)) > 0, ]
  clipped <- suppressWarnings(sf::st_intersection(full, sf::st_geometry(aoi)))
  clipped$area_m2 <- as.numeric(sf::st_area(clipped))
  clipped <- clipped[clipped$area_m2 > 1e-6, ]
  rownames(clipped) <- NULL
  # Polygon intersections give exact AOI area weights, including boundary cells.
  if (abs(sum(clipped$area_m2) / as.numeric(sf::st_area(aoi)) - 1) > 1e-7) {
    stop("Grid intersections do not cover the AOI exactly.")
  }
  list(cells = clipped, cellsize_m = cellsize_m)
}

.grid_points <- function(cells) {
  sf::st_as_sf(sf::st_drop_geometry(cells), coords = c("x_m", "y_m"),
    crs = sf::st_crs(cells), remove = FALSE)
}

.distance_matrix <- function(grid_points, cores) {
  d <- sf::st_distance(grid_points, cores)
  matrix(as.numeric(d), nrow = nrow(grid_points), ncol = nrow(cores))
}

idw_values <- function(distance_m, values, power = 2) {
  distance_m <- as.matrix(distance_m)
  if (ncol(distance_m) != length(values) || !length(values) ||
      any(!is.finite(distance_m) | distance_m < 0) || any(!is.finite(values)) ||
      length(power) != 1 || !is.finite(power) || power <= 0) stop("Invalid IDW inputs.")
  weights <- pmax(distance_m, 1e-9)^(-power)
  answer <- as.numeric(weights %*% values / rowSums(weights))
  exact <- which(rowSums(distance_m < 1e-8) > 0)
  for (i in exact) answer[i] <- mean(values[distance_m[i, ] < 1e-8])
  answer
}

.map_theme <- function() {
  ggplot2::theme_void(base_size = 12) + ggplot2::theme(
    legend.position = "right", legend.title = ggplot2::element_text(size = 11),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(4, 4, 4, 4), strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 12))
}

.stock_map <- function(cells, aoi, value, points = NULL, limits = NULL, panel = NULL,
  difference = FALSE, unit = "kg C/m²") {
  cells$.value <- cells[[value]]
  p <- ggplot2::ggplot(cells) + ggplot2::geom_sf(ggplot2::aes(fill = .value), colour = NA) +
    ggplot2::geom_sf(data = aoi, fill = NA, colour = "grey35", linewidth = .3)
  if (!is.null(points)) p <- p + ggplot2::geom_sf(data = points, shape = 21, fill = "white", size = 2.5)
  if (difference) {
    bound <- max(abs(cells$.value), na.rm = TRUE)
    if (!is.finite(bound) || bound == 0) bound <- 1
    p <- p + ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, limits = c(-bound, bound), name = unit)
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(option = "C", limits = limits, name = unit)
  }
  if (!is.null(panel)) { cells$.panel <- panel; p$data <- cells; p <- p + ggplot2::facet_wrap(~.panel) }
  p + .map_theme()
}

.write_raster <- function(cells, value, path, cellsize_m = 500, units = "kg C/m2") {
  # GeoTIFFs include every cell intersecting the boundary. Exact area totals use
  # clipped polygon areas, not a count of these full raster pixels.
  crs <- sf::st_crs(cells)$wkt
  xmin <- min(cells$x_m) - cellsize_m / 2; xmax <- max(cells$x_m) + cellsize_m / 2
  ymin <- min(cells$y_m) - cellsize_m / 2; ymax <- max(cells$y_m) + cellsize_m / 2
  r <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    resolution = cellsize_m, crs = crs)
  v <- rep(NA_real_, terra::ncell(r))
  v[terra::cellFromXY(r, as.matrix(sf::st_drop_geometry(cells)[c("x_m", "y_m")]))] <- cells[[value]]
  terra::values(r) <- v; names(r) <- value; terra::units(r) <- units
  terra::writeRaster(r, path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  invisible(r)
}

run_area_analysis <- function(config, core_path = "outputs/part1/core_stocks.csv", out = "outputs/part2/area", cellsize_m = 500) {
  x <- .read_spatial_inputs(config, core_path)
  grid <- .area_grid(x$aoi, cellsize_m); cells <- grid$cells
  d <- .distance_matrix(.grid_points(cells), x$points)
  nearest <- max.col(-d, ties.method = "first")
  cells$nearest_core_id <- x$cores$core_id[nearest]
  cells$illustrative_ecosystem <- x$cores$ecosystem[nearest]
  cells$group <- x$cores$group[nearest]
  # One core per ecosystem cannot estimate within-ecosystem sampling variation.
  # A group mean is applied uniformly to each hypothetical forest/wetland area.
  group_mean <- tapply(x$cores$stock_0_30_kg_m2, x$cores$group, mean)
  cells$stock_0_30_kg_m2 <- unname(group_mean[cells$group])
  totals <- do.call(rbind, lapply(split(sf::st_drop_geometry(cells), cells$group), function(z) {
    mean_stock <- unname(group_mean[z$group[1]])
    data.frame(group = z$group[1], area_km2 = sum(z$area_m2) / 1e6,
      area_ha = sum(z$area_m2) / 1e4, n_cores = sum(x$cores$group == z$group[1]),
      mean_stock_0_30_kg_m2 = mean_stock,
      total_stock_0_30_t_C = sum(z$area_m2) * mean_stock / 1000,
      area_basis = paste0("Exact AOI intersections of hypothetical ", cellsize_m, " m nearest-core zones"),
      data_origin = "ILLUSTRATIVE ecosystems and representative group means")
  }))
  rownames(totals) <- NULL
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(totals, file.path(out, "area_totals.csv"), row.names = FALSE)
  utils::write.csv(sf::st_drop_geometry(cells), file.path(out, "grid_values.csv"), row.names = FALSE)
  sf::st_write(cells, file.path(out, "illustrative_ecosystem_cells.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  saveRDS(list(cells = cells, aoi = x$aoi, points = x$points, cores = x$cores,
    cellsize_m = cellsize_m, area_totals = totals), file.path(out, "area_analysis.rds"))
  p <- ggplot2::ggplot() + ggplot2::geom_sf(data = x$aoi, fill = "#f3f3f3", colour = "grey35") +
    ggplot2::geom_sf(data = x$points, ggplot2::aes(colour = ecosystem), size = 3) +
    ggplot2::labs(colour = NULL) + .map_theme()
  ggplot2::ggsave(file.path(out, "study_area.png"), p, width = 7, height = 6, dpi = 180)
  p <- ggplot2::ggplot(cells) + ggplot2::geom_sf(ggplot2::aes(fill = illustrative_ecosystem), colour = NA) +
    ggplot2::geom_sf(data = x$aoi, fill = NA, colour = "grey35", linewidth = .3) +
    ggplot2::geom_sf(data = x$points, shape = 21, fill = "white", size = 2.5) +
    ggplot2::labs(fill = NULL) + .map_theme()
  ggplot2::ggsave(file.path(out, "ecosystem_map.png"), p, width = 8, height = 6, dpi = 180)
  p <- .stock_map(cells, x$aoi, "stock_0_30_kg_m2", x$points)
  ggplot2::ggsave(file.path(out, "group_stock_map.png"), p, width = 7, height = 6, dpi = 180)
  .write_raster(cells, "stock_0_30_kg_m2", file.path(out, "group_stock_0_30.tif"), cellsize_m)
  .spatial_metadata(out, "Hypothetical nearest-core ecosystem areas; group mean upscaling",
    "Total t C = sum(cell area m2 * group mean kg C/m2)/1000",
    "sf; terra; ggplot2; base R",
    "Ecosystems are invented nearest-core zones, not mapped land cover. Four illustrative cores represent whole groups by assumption. Common depth 0–30 cm. Clipped cell areas exactly cover AOI; GeoTIFF boundary pixels are whole cells. No uncertainty or inventory accuracy is inferred.")
  message("Illustrative areas saved to ", out)
  invisible(list(cells = cells, area_totals = totals, aoi = x$aoi))
}

run_spatial_analysis <- function(config, area_path = "outputs/part2/area/area_analysis.rds", out = "outputs/part2/spatial", seed = 20261026) {
  .spatial_packages()
  if (!isTRUE(config$illustrative)) stop("This simulation requires illustrative = TRUE.")
  if (!file.exists(area_path)) stop("Run scripts/03_area.R first.")
  a <- readRDS(area_path); cells <- a$cells
  d <- .distance_matrix(.grid_points(cells), a$points)
  cells$idw_stock_0_30_kg_m2 <- idw_values(d, a$cores$stock_0_30_kg_m2)
  cells$nearest_core_distance_m <- apply(d, 1, min)
  # Entirely synthetic predictors: smooth arbitrary surfaces, not surveyed terrain.
  xn <- (cells$x_m - min(cells$x_m)) / diff(range(cells$x_m))
  yn <- (cells$y_m - min(cells$y_m)) / diff(range(cells$y_m))
  cells$simulated_wetness <- pmax(.01, pmin(.99, .5 + .3 * sin(2*pi*xn) * cos(pi*yn) + .12*yn))
  cells$simulated_elevation_m <- 22 + 14*xn + 5*yn + 3*sin(2*pi*yn)
  set.seed(seed)
  n_train <- min(160L, floor(nrow(cells) / 3))
  sampled <- sample(seq_len(nrow(cells)), n_train)
  simulated <- sf::st_drop_geometry(cells[sampled, c("cell_id", "x_m", "y_m", "simulated_wetness", "simulated_elevation_m")])
  # Anchor only the mean to the revised example. These are additional SIMULATED
  # responses, never extra field samples and never external validation.
  anchor <- mean(a$cores$stock_0_30_kg_m2)
  elevation_mean <- mean(cells$simulated_elevation_m)
  simulated$true_simulated_stock_kg_m2 <- anchor + 5*(simulated$simulated_wetness - .5) -
    .12*(simulated$simulated_elevation_m - elevation_mean)
  simulated$simulated_stock_kg_m2 <- pmax(.1, simulated$true_simulated_stock_kg_m2 + stats::rnorm(n_train, 0, 1.2))
  # A single east/west division prevents random mixing of neighbouring cells.
  cutoff <- stats::quantile(simulated$x_m, .70, names = FALSE)
  simulated$split <- ifelse(simulated$x_m > cutoff, "Spatial holdout", "Training")
  simulated$record_origin <- "SIMULATED; not field observations"
  train <- simulated[simulated$split == "Training", ]; test <- simulated[simulated$split == "Spatial holdout", ]
  fit_train <- stats::lm(simulated_stock_kg_m2 ~ simulated_wetness + simulated_elevation_m, data = train)
  test$predicted_stock_kg_m2 <- as.numeric(stats::predict(fit_train, newdata = test))
  test$residual_kg_m2 <- test$simulated_stock_kg_m2 - test$predicted_stock_kg_m2
  metrics <- data.frame(validation = "Spatial holdout of SIMULATED observations",
    n_training = nrow(train), n_holdout = nrow(test),
    rmse_kg_m2 = sqrt(mean(test$residual_kg_m2^2)), mae_kg_m2 = mean(abs(test$residual_kg_m2)),
    r_squared = if (stats::sd(test$simulated_stock_kg_m2) > 0)
      1-sum(test$residual_kg_m2^2)/sum((test$simulated_stock_kg_m2-mean(test$simulated_stock_kg_m2))^2) else NA_real_,
    applicability = "Simulation performance only; does not validate real community stocks")
  # Refit all simulation records only after holdout evaluation; distinguish this
  # final demonstration map from the held-out model used for the error metrics.
  fit_all <- stats::lm(simulated_stock_kg_m2 ~ simulated_wetness + simulated_elevation_m, data = simulated)
  cells$model_stock_0_30_kg_m2 <- as.numeric(stats::predict(fit_all, newdata = sf::st_drop_geometry(cells)))
  coefficients <- data.frame(term = names(stats::coef(fit_all)), estimate = unname(stats::coef(fit_all)))
  totals <- data.frame(method = c("Group means", "IDW", "Simulated predictor model"),
    total_stock_0_30_t_C = c(sum(cells$area_m2*cells$stock_0_30_kg_m2),
      sum(cells$area_m2*cells$idw_stock_0_30_kg_m2), sum(cells$area_m2*cells$model_stock_0_30_kg_m2))/1000,
    scope = "Illustrative; exact AOI area weighting; 0–30 cm")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (nm in c("simulation_training", "holdout_predictions", "holdout_metrics", "model_coefficients", "method_totals", "grid_values")) {
    z <- switch(nm, simulation_training = simulated, holdout_predictions = test, holdout_metrics = metrics,
      model_coefficients = coefficients, method_totals = totals, grid_values = sf::st_drop_geometry(cells))
    utils::write.csv(z, file.path(out, paste0(nm, ".csv")), row.names = FALSE)
  }
  .write_raster(cells, "idw_stock_0_30_kg_m2", file.path(out, "idw_stock_0_30.tif"), a$cellsize_m)
  .write_raster(cells, "nearest_core_distance_m", file.path(out, "nearest_core_distance_m.tif"), a$cellsize_m, "m")
  .write_raster(cells, "model_stock_0_30_kg_m2", file.path(out, "model_stock_0_30.tif"), a$cellsize_m)
  .write_raster(cells, "simulated_wetness", file.path(out, "simulated_wetness.tif"), a$cellsize_m, "1")
  .write_raster(cells, "simulated_elevation_m", file.path(out, "simulated_elevation_m.tif"), a$cellsize_m, "m")
  ggplot2::ggsave(file.path(out, "idw_stock_map.png"), .stock_map(cells, a$aoi, "idw_stock_0_30_kg_m2", a$points), width = 7, height = 6, dpi = 180)
  ggplot2::ggsave(file.path(out, "distance_map.png"), .stock_map(cells, a$aoi, "nearest_core_distance_m", a$points, unit = "m"), width = 7, height = 6, dpi = 180)
  ggplot2::ggsave(file.path(out, "model_stock_map.png"), .stock_map(cells, a$aoi, "model_stock_0_30_kg_m2", a$points), width = 7, height = 6, dpi = 180)
  plot_limits <- range(c(test$simulated_stock_kg_m2, test$predicted_stock_kg_m2))
  p <- ggplot2::ggplot(test, ggplot2::aes(simulated_stock_kg_m2, predicted_stock_kg_m2)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey55", linetype = 2) +
    ggplot2::geom_point(colour = "#317873", size = 2) +
    ggplot2::coord_equal(xlim = plot_limits, ylim = plot_limits) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::labs(x = "Simulated (kg C/m²)", y = "Predicted (kg C/m²)")
  ggplot2::ggsave(file.path(out, "observed_predicted.png"), p, width = 6, height = 5, dpi = 180)
  saveRDS(list(cells = cells, aoi = a$aoi, points = a$points, cores = a$cores,
    cellsize_m = a$cellsize_m, simulation = simulated, holdout = test, metrics = metrics,
    training_model = fit_train, final_model = fit_all), file.path(out, "spatial_analysis.rds"))
  .spatial_metadata(out,
    c("IDW", "Distance to nearest core", "Synthetic two-predictor regression", "Spatial holdout"),
    c("stock(s) = sum(d_i^-2 * stock_i)/sum(d_i^-2); exact input value when d_i = 0",
      "distance(s) = min_i d(s, core_i)",
      "stats::lm(simulated_stock_kg_m2 ~ simulated_wetness + simulated_elevation_m)",
      "RMSE = sqrt(mean((simulated-predicted)^2)); MAE = mean(abs(simulated-predicted))"),
    c("base R; sf::st_distance; terra; ggplot2", "base R; sf::st_distance; terra; ggplot2", "stats::lm; terra; ggplot2", "stats::predict; base R; ggplot2"),
    c("Four illustrative core stocks; distance power 2; no variogram or kriging uncertainty. Common 0–30 cm depth.",
      "Proximity only; not a variance, confidence interval or map accuracy estimate.",
      "Separate simulated predictors and responses; their mean uses the four-core mean. These are not additional ground data. Final model refit to all simulation records after evaluation.",
      "Eastern approximately 30% held out by location. Simulation performance only; not independent validation of field data or map accuracy."))
  message("Illustrative spatial outputs saved to ", out)
  invisible(list(cells = cells, metrics = metrics, simulation = simulated, holdout = test, totals = totals))
}

normal_update <- function(prior_mean, prior_sd, new_mean, new_sd) {
  if (any(!is.finite(c(prior_mean, prior_sd, new_mean, new_sd))) ||
      any(prior_sd <= 0) || any(new_sd <= 0)) stop("Normal update requires finite means and positive SDs.")
  lengths <- c(length(prior_mean), length(prior_sd), length(new_mean), length(new_sd))
  if (any(!lengths) || any(!lengths %in% c(1L, max(lengths)))) stop("Normal update lengths must agree or be scalar.")
  variance <- 1 / (1/prior_sd^2 + 1/new_sd^2)
  mean <- variance * (prior_mean/prior_sd^2 + new_mean/new_sd^2)
  # Scalar variance is deliberately recycled across all matching mean values.
  data.frame(mean = mean, variance = rep_len(variance, length(mean)),
    sd = rep_len(sqrt(variance), length(mean)))
}

.save_map_panels <- function(plots, path, width = 13, height = 5) {
  grDevices::png(path, width = width, height = height, units = "in", res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(plots))))
  for (i in seq_along(plots)) print(plots[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
}

run_bayesian_update <- function(config, spatial_path = "outputs/part2/spatial/spatial_analysis.rds", out = "outputs/part3", prior_sd = 5, new_sd = 3) {
  .spatial_packages()
  if (!isTRUE(config$illustrative)) stop("This example requires illustrative = TRUE.")
  if (!file.exists(spatial_path)) stop("Run scripts/04_spatial.R first.")
  a <- readRDS(spatial_path); cells <- a$cells
  xn <- (cells$x_m - min(cells$x_m)) / diff(range(cells$x_m))
  # Arbitrary prior map, created without using workbook stocks or simulation
  # residuals. Both uncertainty levels are teaching assumptions, not estimates.
  cells$prior_mean_kg_m2 <- 12 + 2*(xn - .5)
  cells$prior_sd_kg_m2 <- prior_sd
  cells$new_mean_kg_m2 <- cells$idw_stock_0_30_kg_m2
  cells$new_sd_kg_m2 <- new_sd
  post <- normal_update(cells$prior_mean_kg_m2, cells$prior_sd_kg_m2,
    cells$new_mean_kg_m2, cells$new_sd_kg_m2)
  cells$posterior_mean_kg_m2 <- post$mean; cells$posterior_variance <- post$variance
  cells$posterior_sd_kg_m2 <- post$sd
  cells$posterior_lower95_kg_m2 <- post$mean - stats::qnorm(.975)*post$sd
  cells$posterior_upper95_kg_m2 <- post$mean + stats::qnorm(.975)*post$sd
  cells$change_kg_m2 <- cells$posterior_mean_kg_m2 - cells$prior_mean_kg_m2
  cells$new_weight <- cells$prior_sd_kg_m2^2/(cells$prior_sd_kg_m2^2 + cells$new_sd_kg_m2^2)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(sf::st_drop_geometry(cells), file.path(out, "bayesian_grid.csv"), row.names = FALSE)
  summary <- data.frame(stage = c("Prior", "New estimate", "Updated"),
    area_weighted_mean_kg_m2 = c(weighted.mean(cells$prior_mean_kg_m2, cells$area_m2),
      weighted.mean(cells$new_mean_kg_m2, cells$area_m2), weighted.mean(cells$posterior_mean_kg_m2, cells$area_m2)),
    total_stock_0_30_t_C = c(sum(cells$prior_mean_kg_m2*cells$area_m2),
      sum(cells$new_mean_kg_m2*cells$area_m2), sum(cells$posterior_mean_kg_m2*cells$area_m2))/1000,
    cell_sd_kg_m2 = c(prior_sd, new_sd, sqrt(1/(1/prior_sd^2 + 1/new_sd^2))),
    interpretation = "Illustrative cell-level uncertainty; not an uncertainty interval for the AOI total")
  utils::write.csv(summary, file.path(out, "bayesian_summary.csv"), row.names = FALSE)
  for (nm in c("prior_mean_kg_m2", "prior_sd_kg_m2", "new_mean_kg_m2", "new_sd_kg_m2",
    "posterior_mean_kg_m2", "posterior_sd_kg_m2", "posterior_lower95_kg_m2", "posterior_upper95_kg_m2", "change_kg_m2")) {
    .write_raster(cells, nm, file.path(out, paste0(nm, ".tif")), a$cellsize_m)
  }
  limits <- range(c(cells$prior_mean_kg_m2, cells$posterior_mean_kg_m2))
  plots <- list(.stock_map(cells, a$aoi, "prior_mean_kg_m2", limits = limits, panel = "Prior"),
    .stock_map(cells, a$aoi, "posterior_mean_kg_m2", limits = limits, panel = "Updated"),
    .stock_map(cells, a$aoi, "change_kg_m2", panel = "Change", difference = TRUE))
  .save_map_panels(plots, file.path(out, "bayesian_maps.png"))
  p <- ggplot2::ggplot(summary, ggplot2::aes(factor(stage, levels = c("Prior", "New estimate", "Updated")), cell_sd_kg_m2)) +
    ggplot2::geom_col(fill = "#317873", width = .6) + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = NULL, y = "SD (kg C/m²)")
  ggplot2::ggsave(file.path(out, "uncertainty.png"), p, width = 6, height = 4, dpi = 180)
  saveRDS(list(cells = cells, summary = summary), file.path(out, "bayesian_analysis.rds"))
  .spatial_metadata(out, "Normal–normal illustrative update",
    "V_post = 1/(1/V_prior + 1/V_new); mu_post = V_post*(mu_prior/V_prior + mu_new/V_new); SD_post = sqrt(V_post)",
    "base R; stats::qnorm; sf; terra; ggplot2",
    paste0("Prior mean = 12 + 2*(scaled east coordinate - 0.5) kg C/m2, independent of workbook inputs. Prior SD = ", prior_sd,
      "; new estimate = four-core IDW; assumed new-estimate SD = ", new_sd,
      ". Errors treated as independent Gaussian cell-level errors. SDs are assumptions, not between-core SDs or fitted standard errors. Shared spatial errors and model bias are omitted. Conditional intervals do not validate real stocks; no AOI total uncertainty is calculated. Normal support may include negative stock in other examples."))
  message("Illustrative Bayesian outputs saved to ", out)
  invisible(list(cells = cells, summary = summary))
}
