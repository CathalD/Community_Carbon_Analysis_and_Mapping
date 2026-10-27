# Spatial methods with offline geographic context.
source("R/cartography.R")
# All stock comparisons use soil organic carbon at 0–30 cm, in kg C/m2.

.spatial_packages <- function() {
  for (p in c("sf", "terra", "ggplot2", "gstat")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, " before this step.")
  }
}

.spatial_metadata <- function(out, method, formula, packages, assumptions) {
  utils::write.csv(data.frame(method = method, formula = formula,
    packages = packages, assumptions = assumptions),
    file.path(out, "method_metadata.csv"), row.names = FALSE)
}

.read_spatial_inputs <- function(config, core_path = "outputs/part1/core_stocks.csv", excluded_path = NULL) {
  .spatial_packages()
  if (!file.exists(core_path)) stop("Calculate the field stocks first; core_stocks.csv is missing.")
  if (is.null(config$aoi_path) || !file.exists(config$aoi_path)) {
    stop("AOI file not found. Save your boundary as data/aoi/my_aoi.geojson and check aoi_path in your settings.")
  }
  cores <- utils::read.csv(core_path, stringsAsFactors = FALSE)
  needed <- c("core_id", "latitude", "longitude", "stock_0_30_kg_m2")
  if (!all(needed %in% names(cores))) stop("Core stock table lacks spatial input columns.")
  if (anyDuplicated(cores$core_id)) stop("Core IDs must be unique.")
  for (v in c("latitude", "longitude", "stock_0_30_kg_m2")) {
    cores[[v]] <- suppressWarnings(as.numeric(cores[[v]]))
  }
  for (v in c("ecosystem", "group")) if (!v %in% names(cores)) cores[[v]] <- "Not recorded"
  aoi <- sf::st_read(config$aoi_path, quiet = TRUE)
  if (!nrow(aoi) || any(sf::st_is_empty(aoi))) stop("The AOI boundary is empty. Export a closed polygon boundary.")
  if (!all(as.character(sf::st_geometry_type(aoi)) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("The AOI must be a polygon boundary, not points or lines. Export it as a polygon GeoJSON.")
  }
  if (is.na(sf::st_crs(aoi))) stop("The AOI needs a coordinate reference system. Export WGS84 GeoJSON from your GIS.")
  aoi <- sf::st_transform(sf::st_make_valid(aoi), 4326)
  bb <- sf::st_bbox(aoi)
  if (any(!is.finite(bb))) stop("The AOI has invalid coordinates. Check its coordinate reference system.")
  if (bb["xmax"] - bb["xmin"] > 180) stop("AOI crosses the antimeridian or has wrapped longitudes. Split and normalise it in a GIS first.")
  if (bb["xmin"] < -170 || bb["xmax"] > -40 || bb["ymin"] < 35 || bb["ymax"] > 88)
    warning("AOI extends beyond the Canadian context cache. Check longitude signs and CRS; context may be incomplete.")
  local_crs <- paste0("+proj=laea +lat_0=", (bb["ymin"] + bb["ymax"]) / 2,
    " +lon_0=", (bb["xmin"] + bb["xmax"]) / 2, " +datum=WGS84 +units=m +no_defs")
  aoi <- sf::st_transform(aoi, local_crs)
  if (any(!as.character(sf::st_geometry_type(aoi)) %in% c("POLYGON", "MULTIPOLYGON"))) {
    aoi <- tryCatch(suppressWarnings(sf::st_collection_extract(aoi, "POLYGON")),
      error = function(e) stop("The AOI does not enclose a valid polygon. Check its boundary in your GIS."))
  }
  aoi <- sf::st_sf(area_id = "Study area", geometry = sf::st_union(sf::st_make_valid(aoi)))
  area <- as.numeric(sf::st_area(aoi))
  if (length(area) != 1 || !is.finite(area) || area <= 0) stop("The AOI must enclose an area greater than zero.")
  projected_bounds <- sf::st_bbox(aoi)
  if (max(projected_bounds["xmax"]-projected_bounds["xmin"], projected_bounds["ymax"]-projected_bounds["ymin"]) > 1e6)
    warning("AOI spans over 1,000 km: LAEA preserves area but distance distortion matters. Treat this IDW as a technical demonstration and split into regional analyses for scientific use.")
  reason <- rep("", nrow(cores))
  add_reason <- function(which_rows, text) {
    reason[which_rows] <<- ifelse(nzchar(reason[which_rows]), paste(reason[which_rows], text, sep = "; "), text)
  }
  valid_coordinates <- is.finite(cores$latitude) & is.finite(cores$longitude) &
    abs(cores$latitude) <= 90 & abs(cores$longitude) <= 180
  add_reason(!valid_coordinates, "Missing or invalid WGS84 latitude/longitude")
  add_reason(!is.finite(cores$stock_0_30_kg_m2), "Incomplete 0–30 cm stock")
  located <- which(valid_coordinates)
  if (length(located)) {
    located_points <- sf::st_as_sf(cores[located, ], coords = c("longitude", "latitude"), crs = 4326)
    inside <- lengths(sf::st_intersects(sf::st_transform(located_points, local_crs), aoi)) > 0
    add_reason(located[!inside], "Outside the AOI boundary")
  }
  excluded <- cores[nzchar(reason), , drop = FALSE]
  excluded$reason <- reason[nzchar(reason)]
  if (!is.null(excluded_path)) utils::write.csv(excluded, excluded_path, row.names = FALSE, na = "")
  if (nrow(excluded)) message("Excluded ", nrow(excluded), " core(s) from mapping: ",
    paste(excluded$core_id, collapse = ", "), ". Review excluded_cores.csv; field results are retained.")
  cores <- cores[!nzchar(reason), , drop = FALSE]
  if (nrow(cores) < 2) {
    stop("Mapping needs at least two cores at different locations inside the AOI, with complete 0–30 cm stocks. Check coordinates, section depths and your boundary; field outputs are still available.")
  }
  if (anyDuplicated(cores[c("latitude", "longitude")])) {
    stop("Two mapped cores share the same coordinates. For this simple map, keep one representative core per location or combine repeated cores before mapping.")
  }
  points <- sf::st_as_sf(cores, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  points <- sf::st_transform(points, local_crs)
  list(aoi = aoi, points = points, cores = cores, excluded = excluded)
}

.area_grid <- function(aoi, cellsize_m = 500) {
  bb <- sf::st_bbox(aoi)
  span <- c(unname(bb["xmax"] - bb["xmin"]), unname(bb["ymax"] - bb["ymin"]))
  if (any(!is.finite(span)) || any(span <= 0)) stop("The AOI must have a nonempty polygon boundary.")
  if (is.null(cellsize_m)) {
    # Bounding-box area also protects against widely separated AOI polygons.
    cellsize_m <- max(sqrt(prod(span) / 2500), max(span) / 2500, .001)
    cellsize_m <- signif(cellsize_m, 2)
  }
  if (length(cellsize_m) != 1 || !is.finite(cellsize_m) || cellsize_m <= 0) stop("Invalid grid resolution.")
  candidates <- prod(pmax(1, ceiling(span / cellsize_m)))
  if (!is.finite(candidates) || candidates > 10000) {
    stop("This grid would have more than 10,000 cells. Increase cellsize_m, or use NULL for an automatic grid.")
  }
  grid <- sf::st_make_grid(aoi, cellsize = cellsize_m, square = TRUE)
  full <- sf::st_sf(cell_id = seq_along(grid), geometry = grid)
  centres <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(full)))
  full$x_m <- centres[, 1]; full$y_m <- centres[, 2]
  full <- full[lengths(sf::st_intersects(full, aoi)) > 0, ]
  clipped <- suppressWarnings(sf::st_intersection(full, sf::st_geometry(aoi)))
  clipped$area_m2 <- as.numeric(sf::st_area(clipped))
  clipped <- clipped[is.finite(clipped$area_m2) & clipped$area_m2 > 0, ]
  rownames(clipped) <- NULL
  # Polygon intersections give exact AOI area weights, including boundary cells.
  if (abs(sum(clipped$area_m2) / as.numeric(sf::st_area(aoi)) - 1) > 1e-7) {
    stop("Grid intersections do not cover the AOI exactly.")
  }
  list(cells = clipped, cellsize_m = cellsize_m, candidate_cells = candidates)
}

.grid_points <- function(cells) {
  sf::st_as_sf(sf::st_drop_geometry(cells), coords = c("x_m", "y_m"),
    crs = sf::st_crs(cells), remove = FALSE)
}

.distance_matrix <- function(grid_points, cores) {
  if (as.double(nrow(grid_points)) * nrow(cores) > 2e6)
    stop("Distance matrix exceeds two million pairs. Increase grid resolution or subset observations.")
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
  if (is.null(limits) && !difference && unit == "kg C/m²" && !is.null(points))
    limits <- range(c(cells[[value]], points$stock_0_30_kg_m2), na.rm = TRUE)
  cells$.value <- cells[[value]]
  p <- .context_base(aoi) + ggplot2::geom_sf(data = cells, ggplot2::aes(fill = .value), colour = NA) +
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
  .context_finish(p, aoi)
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

run_simple_mapping <- function(config, core_path = "outputs/part1/core_stocks.csv",
  out = "outputs/my_data/maps", cellsize_m = NULL) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  paths <- list(location_map = file.path(out, "location_map.png"), sample_map = file.path(out, "sample_map.png"),
    idw_map = file.path(out, "idw_stock_map.png"),
    idw_raster = file.path(out, "idw_stock_0_30.tif"),
    grid_values = file.path(out, "grid_values.csv"),
    metadata = file.path(out, "mapping_metadata.csv"),
    excluded_cores = file.path(out, "excluded_cores.csv"))
  # A failed rerun must not leave yesterday's maps beside today's field data.
  # Remove only this function's known outputs; retain all other folder contents.
  previous_outputs <- unlist(paths, use.names = FALSE)
  previous_outputs <- previous_outputs[file.exists(previous_outputs)]
  if (length(previous_outputs) && unlink(previous_outputs) != 0) {
    stop("Could not clear an earlier mapping output. Close that file and run the mapping step again.")
  }
  x <- .read_spatial_inputs(config, core_path, paths$excluded_cores)
  ggplot2::ggsave(paths$location_map, regional_location_map(x$aoi), width=8, height=6, dpi=180)
  grid <- .area_grid(x$aoi, cellsize_m)
  cells <- grid$cells
  distance <- .distance_matrix(.grid_points(cells), x$points)
  cells$idw_stock_0_30_kg_m2 <- .gstat_idw(.grid_points(cells), x$points)
  cells$nearest_core_distance_m <- apply(distance, 1, min)
  p <- .context_base(x$aoi) +
    ggplot2::geom_sf(data = x$aoi, fill = "#f3f3f3", colour = "grey35", linewidth = .3) +
    ggplot2::geom_sf(data = x$points, ggplot2::aes(fill = stock_0_30_kg_m2),
      shape = 21, colour = "grey20", size = 3) +
    ggplot2::scale_fill_viridis_c(option = "C", name = "kg C/m²") + .map_theme()
  ggplot2::ggsave(paths$sample_map, .context_finish(p, x$aoi), width = 7, height = 6, dpi = 180)
  ggplot2::ggsave(paths$idw_map,
    .stock_map(cells, x$aoi, "idw_stock_0_30_kg_m2", x$points), width = 7, height = 6, dpi = 180)
  .write_raster(cells, "idw_stock_0_30_kg_m2", paths$idw_raster, grid$cellsize_m)
  utils::write.csv(sf::st_drop_geometry(cells), paths$grid_values, row.names = FALSE)
  bb <- sf::st_bbox(sf::st_transform(x$aoi, 4326))
  metadata <- data.frame(method = "Inverse distance weighting (IDW), power 2",
    packages = "gstat::idw; sf; terra; ggplot2",
    formula = "stock(s) = sum(d_i^-2 * stock_i) / sum(d_i^-2)",
    stock_depth_cm = "0–30", unit = "kg C/m2", n_cores = nrow(x$cores),
    excluded_cores = nrow(x$excluded), grid_cellsize_m = grid$cellsize_m,
    candidate_grid_cells = grid$candidate_cells, mapped_grid_cells = nrow(cells),
    aoi_area_km2 = as.numeric(sf::st_area(x$aoi)) / 1e6,
    longitude_min = unname(bb["xmin"]), longitude_max = unname(bb["xmax"]),
    latitude_min = unname(bb["ymin"]), latitude_max = unname(bb["ymax"]),
    source_coordinates = "WGS84 decimal-degree latitude and longitude",
    map_projection = sf::st_crs(x$aoi)$input,
    data_origin = if (is.null(config$data_origin)) "User-supplied data" else config$data_origin,
    interpretation = paste("Nearby cores influence the estimate more strongly.",
      "IDW shows a smooth estimate; it does not measure map accuracy or uncertainty.",
      "No ecosystem areas, extra observations, prior maps or AOI carbon inventory are invented.",
      "GeoTIFF boundary pixels are whole cells; larger AOIs use coarser grids."))
  utils::write.csv(metadata, paths$metadata, row.names = FALSE)
  message("Mapped ", nrow(x$cores), " core(s) at ", signif(grid$cellsize_m, 3),
    " m resolution. Outputs: ", out)
  invisible(list(cells = cells, aoi = x$aoi, points = x$points, cores = x$cores,
    excluded = x$excluded, cellsize_m = grid$cellsize_m, metadata = metadata, paths = paths))
}

run_area_analysis <- function(config, core_path = "outputs/part1/core_stocks.csv", out = "outputs/part2/area", cellsize_m = 500) {
  if (!isTRUE(config$illustrative)) stop("This ecosystem-area example requires illustrative = TRUE. Use run_simple_mapping() for your own measurements.")
  x <- .read_spatial_inputs(config, core_path)
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(file.path(out,"location_map.png"), regional_location_map(x$aoi), width=8,height=6,dpi=180)
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
  p <- .context_base(x$aoi) + ggplot2::geom_sf(data = x$aoi, fill = "#f3f3f3", colour = "grey35") +
    ggplot2::geom_sf(data = x$points, ggplot2::aes(colour = ecosystem), size = 3) +
    ggplot2::labs(colour = NULL) + .map_theme()
  ggplot2::ggsave(file.path(out, "study_area.png"), .context_finish(p, x$aoi), width = 7, height = 6, dpi = 180)
  p <- .context_base(x$aoi) + ggplot2::geom_sf(data = cells, ggplot2::aes(fill = illustrative_ecosystem), colour = NA) +
    ggplot2::geom_sf(data = x$aoi, fill = NA, colour = "grey35", linewidth = .3) +
    ggplot2::geom_sf(data = x$points, shape = 21, fill = "white", size = 2.5) +
    ggplot2::labs(fill = NULL) + .map_theme()
  ggplot2::ggsave(file.path(out, "ecosystem_map.png"), .context_finish(p, x$aoi), width = 8, height = 6, dpi = 180)
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
  cells$idw_stock_0_30_kg_m2 <- .gstat_idw(.grid_points(cells), a$points)
  cells$nearest_core_distance_m <- apply(d, 1, min)
  # Entirely synthetic predictors: smooth arbitrary surfaces, not surveyed terrain.
  xn <- .scale_coordinate(cells$x_m)
  yn <- .scale_coordinate(cells$y_m)
  cells$simulated_wetness <- pmax(.01, pmin(.99, .5 + .3 * sin(2*pi*xn) * cos(pi*yn) + .12*yn))
  cells$simulated_elevation_m <- 22 + 14*xn + 5*yn + 3*sin(2*pi*yn)
  if (nrow(cells) < 60 || length(unique(cells$x_m)) < 4 || length(unique(cells$y_m)) < 4)
    stop("Landscape simulation needs at least 60 cells and four rows and columns; use a finer grid.")
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
  if (nrow(train) < 10 || nrow(test) < 5) stop("Spatial split has too few training or test locations.")
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
  prior <- 12 + 2*(.scale_coordinate(cells$x_m)-.5)
  prior_comparison <- data.frame(method=c("Group means", "IDW", "Simulated predictor model"),
    prior_total_t_C=sum(prior*cells$area_m2)/1000,
    estimate_total_t_C=totals$total_stock_0_30_t_C,
    difference_t_C=totals$total_stock_0_30_t_C-sum(prior*cells$area_m2)/1000,
    scope="Invented independent prior; same AOI, point-grid approximation and 0-30 cm surface stock; no accuracy ranking")
  write.csv(prior_comparison,file.path(out,"prior_comparison.csv"),row.names=FALSE)
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
  p <- ggplot2::ggplot(test, ggplot2::aes(predicted_stock_kg_m2, simulated_stock_kg_m2)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey55", linetype = 2) +
    ggplot2::geom_point(colour = "#317873", size = 2) +
    ggplot2::coord_equal(xlim = plot_limits, ylim = plot_limits) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::labs(x = "Predicted (kg C/m²)", y = "Simulated observed (kg C/m²)")
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
    c("gstat::idw; sf::st_distance; terra; ggplot2", "base R; sf::st_distance; terra; ggplot2", "stats::lm; terra; ggplot2", "stats::predict; base R; ggplot2"),
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
  xn <- .scale_coordinate(cells$x_m)
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

# Production IDW uses gstat; idw_values remains a transparent arithmetic check.
.gstat_idw <- function(targets, points, power = 2) {
  as.numeric(gstat::idw(stock_0_30_kg_m2 ~ 1, points, targets,
    idp = power, debug.level = 0)$var1.pred)
}
.scale_coordinate <- function(x) {
  span <- diff(range(x))
  if (!is.finite(span) || span == 0) return(rep(.5, length(x)))
  (x - min(x)) / span
}
