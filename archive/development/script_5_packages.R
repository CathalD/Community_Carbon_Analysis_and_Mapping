#==========================================
#----- Script 5: Package-based updates -----
#==========================================
# Run from the workshop project folder.
# The handwritten scripts and outputs are only read, never sourced or overwritten.
# gstat and geoR replace spatial conditioning; mvtnorm replaces scalar conditioning.
# Model assumptions, transformations and covariance specification remain explicit.

library(terra)
library(readr)
library(ggplot2)
terra::projNetwork(FALSE)

required_packages <- c("gstat", "geoR", "mvtnorm", "sp", "sf")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
if (!all(c("cond_mvnorm", "syMatrices", "Tcrossprod") %in% getNamespaceExports("mvtnorm"))) {
  stop("Update mvtnorm: this script needs cond_mvnorm, syMatrices and Tcrossprod.")
}


#==========================================
#----- Settings and inputs ----------------
#==========================================

target_depth_cm <- 100
range_m <- 2000
observation_log_sd <- 0.00
compare_handwritten <- TRUE

stopifnot(target_depth_cm == 100, range_m > 0, observation_log_sd >= 0)

stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
output_folder <- file.path("outputs/Testing_Bayes_Packages", "stock_0_100")
handwritten_folder <- file.path("outputs/Testing_Bayes_Models", "stock_0_100")

input_files <- c(
  field = "outputs/script1/mapping_data.csv",
  boundary = "data/boundary/aoi.geojson",
  prior_mean = "outputs/script4/mcmaster_0_100cm.tif",
  prior_uncertainty = "outputs/script4/mcmaster_uncertainty_0_100cm.tif",
  area_based = "outputs/script2/area_based.gpkg",
  idw = "outputs/script2/idw_0_100cm.tif",
  idw_parameters = "outputs/script2/idw_best_parameters_by_depth.csv"
)

if (!all(file.exists(input_files))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


#==========================================
#----- Prepare inputs independently -------
#==========================================

field <- read_csv(input_files[["field"]], show_col_types = FALSE)
stopifnot(all(c("Core ID", "Plot ID", "Longitude", "Latitude", stock_column) %in% names(field)))
cores <- data.frame(
  core_id = as.character(field[["Core ID"]]),
  plot_id = as.character(field[["Plot ID"]]),
  longitude = field[["Longitude"]], latitude = field[["Latitude"]],
  stock = field[[stock_column]]
)
stopifnot(
  nrow(cores) >= 3, !anyNA(cores), !anyDuplicated(cores$core_id),
  all(nzchar(cores$core_id)), all(nzchar(cores$plot_id)),
  all(is.finite(cores$stock) & cores$stock > 0),
  all(is.finite(cores$longitude) & abs(cores$longitude) <= 180),
  all(is.finite(cores$latitude) & abs(cores$latitude) <= 90),
  length(unique(cores$plot_id)) >= 3
)

aoi <- project(vect(input_files[["boundary"]]), "EPSG:3978")
template <- rast(
  xmin = floor(xmin(aoi) / 1000) * 1000,
  xmax = ceiling(xmax(aoi) / 1000) * 1000,
  ymin = floor(ymin(aoi) / 1000) * 1000,
  ymax = ceiling(ymax(aoi) / 1000) * 1000,
  resolution = 1000, crs = crs(aoi)
)

prior_mean <- rast(input_files[["prior_mean"]])
prior_width <- rast(input_files[["prior_uncertainty"]])
stopifnot(compareGeom(prior_mean, prior_width, stopOnError = FALSE))

# Source uncertainty is Q95-Q05, not SD. This conversion is a working approximation.
native_prior <- c(prior_mean, prior_width / (2 * qnorm(0.95)))
names(native_prior) <- c("prior_mean", "prior_sd")

# Copy paired mean/SD values to grid centres; no block-uncertainty aggregation.
grid_prior <- mask(project(native_prior, template, method = "near"), aoi, touches = FALSE)
usable <- is.finite(grid_prior[[1]]) & grid_prior[[1]] > 0 &
  is.finite(grid_prior[[2]]) & grid_prior[[2]] > 0
grid_prior <- mask(grid_prior, usable, maskvalues = 0)

core_points <- project(
  vect(cores, geom = c("longitude", "latitude"), crs = "EPSG:4326"),
  crs(grid_prior)
)
stopifnot(all(rowSums(relate(core_points, aoi, "intersects")) > 0))
at_cores <- extract(native_prior, project(core_points, crs(native_prior)), ID = FALSE)
stopifnot(all(is.finite(as.matrix(at_cores))), all(as.matrix(at_cores) > 0))
cores$prior_mean <- at_cores$prior_mean
cores$prior_sd <- at_cores$prior_sd
core_xy <- crds(core_points)
if (observation_log_sd == 0 && anyDuplicated(core_xy)) {
  stop("Coincident cores with zero observation error require explicit replicate handling.")
}

area_polygon <- project(vect(input_files[["area_based"]], layer = "area_based_0_100"), crs(template))
area_grid <- mask(
  rasterize(area_polygon, template, field = "mean_stock_kg_m2", touches = FALSE),
  grid_prior[[1]]
)
names(area_grid) <- "field_estimate"
idw_grid <- mask(project(rast(input_files[["idw"]]), template, method = "bilinear"), grid_prior[[1]])
names(idw_grid) <- "field_estimate"

writeRaster(grid_prior, file.path(output_folder, "prepared_prior.tif"), overwrite = TRUE)


#==========================================
#----- Output helpers ---------------------
#==========================================

pixel_table <- function(field_map = NULL) {
  layers <- if (is.null(field_map)) grid_prior else c(grid_prior, field_map)
  d <- as.data.frame(layers, cells = TRUE, xy = TRUE, na.rm = TRUE)
  keep <- is.finite(d$prior_mean) & d$prior_mean > 0 & is.finite(d$prior_sd) & d$prior_sd > 0
  if (!is.null(field_map)) keep <- keep & is.finite(d$field_estimate) & d$field_estimate > 0
  d <- d[keep, , drop = FALSE]
  stopifnot(nrow(d) > 0)
  d
}

log_prior <- function(m, s) {
  v <- log1p((s / m)^2)
  list(mean = log(m) - v / 2, variance = v)
}

stock_results <- function(d, log_mean, log_var) {
  p <- log_prior(d$prior_mean, d$prior_sd)
  stopifnot(all(is.finite(log_mean)), all(is.finite(log_var)))
  if (any(log_var < -1e-8 * pmax(1, p$variance))) stop("Invalid package posterior variance.")
  log_var <- pmax(0, log_var)
  stopifnot(all(log_var <= p$variance + 1e-8))
  d$posterior_log_mean <- log_mean
  d$posterior_log_var <- log_var
  d$posterior_mean <- exp(log_mean + log_var / 2)
  d$posterior_median <- exp(log_mean)
  d$posterior_sd <- sqrt(expm1(log_var) * exp(2 * log_mean + log_var))
  d$lower95 <- exp(log_mean - qnorm(0.975) * sqrt(log_var))
  d$upper95 <- exp(log_mean + qnorm(0.975) * sqrt(log_var))
  d$interval_width95 <- d$upper95 - d$lower95
  d$change_from_prior <- d$posterior_mean - d$prior_mean
  d$log_variance_reduction <- 1 - log_var / p$variance
  stopifnot(all(is.finite(as.matrix(d))))
  d
}

save_analysis <- function(d, method) {
  folder <- file.path(output_folder, method)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  columns <- c("prior_mean", "prior_sd", "posterior_mean", "posterior_median",
               "posterior_sd", "lower95", "upper95", "interval_width95",
               "change_from_prior", "log_variance_reduction")
  maps <- rast(template, nlyrs = length(columns))
  names(maps) <- columns
  vals <- matrix(NA_real_, nrow = ncell(maps), ncol = length(columns))
  vals[d$cell, ] <- as.matrix(d[, columns])
  values(maps) <- vals
  writeRaster(maps, file.path(folder, "posterior_maps.tif"), overwrite = TRUE)
  write_csv(d, file.path(folder, "predictions.csv"))
  plot <- ggplot(d, aes(x, y, fill = posterior_mean)) + geom_raster() +
    scale_fill_viridis_c(name = "kg C/m²") + coord_equal() + theme_minimal() +
    labs(title = paste(method, "package update: 0-100 cm"), x = "Easting (m)", y = "Northing (m)")
  ggsave(file.path(folder, "posterior_mean.png"), plot, width = 8, height = 7, dpi = 300)
}


#==========================================
#----- 1. Raw cores: gstat -----------------
#==========================================

raw_pixels <- pixel_table()
raw_prior <- log_prior(raw_pixels$prior_mean, raw_pixels$prior_sd)
core_prior <- log_prior(cores$prior_mean, cores$prior_sd)

# Standardise the latent process to variance 1. This preserves the original
# varying prior variances after back-transformation and makes krige0's
# constant target-variance assumption appropriate.
projection <- sp::CRS(SRS_string = crs(template))
training <- sp::SpatialPointsDataFrame(
  core_xy,
  data = data.frame(
    residual = (log(cores$stock) - core_prior$mean) / sqrt(core_prior$variance),
    noise_var = observation_log_sd^2 / core_prior$variance
  ),
  proj4string = projection
)
targets <- sp::SpatialPointsDataFrame(
  as.matrix(raw_pixels[, c("x", "y")]),
  data = data.frame(noise_var = rep(0, nrow(raw_pixels))),
  proj4string = projection
)

# Explicit prior covariance specification; gstat performs the conditioning.
# Noise is added to the observation diagonal, not to cross-covariances or targets.
raw_covariance <- function(x, y = x) {
  C <- exp(-sp::spDists(sp::coordinates(x), sp::coordinates(y)) / range_m)
  if (identical(x, y)) diag(C) <- diag(C) + x$noise_var
  C
}

raw_fit <- gstat::krige0(
  residual ~ 1, data = training, newdata = targets,
  model = raw_covariance, beta = 0,
  computeVar = TRUE, fullCovariance = FALSE
)

raw_results <- stock_results(
  raw_pixels,
  raw_prior$mean + sqrt(raw_prior$variance) * as.vector(raw_fit$pred),
  raw_prior$variance * as.vector(raw_fit$var)
)
save_analysis(raw_results, "raw_cores")


#==========================================
#----- 1B. Raw cores: geoR -----------------
#==========================================

# Use the same standardised residual process as gstat: known mean 0,
# variance 1, exponential correlation and fixed range. No range/noise fitting.
# geoR accepts a constant relative nugget, whereas the standardised observation
# noise would vary by core if a common nonzero log-SD were introduced.
geor_noise <- training$noise_var
if (any(abs(geor_noise - geor_noise[1]) > 1e-12)) {
  stop("This geoR comparison needs constant standardised observation variance. ",
       "The current zero-error model qualifies; nonzero heterogeneous noise ",
       "requires a different geoR model. Do not average the noise variances.")
}

geor_fit <- geoR::krige.bayes(
  coords = core_xy,
  data = training$residual,
  locations = as.matrix(raw_pixels[, c("x", "y")]),
  model = geoR::model.control(cov.model = "exponential", lambda = 1),
  prior = geoR::prior.control(
    beta.prior = "fixed", beta = 0,
    sigmasq.prior = "fixed", sigmasq = 1,
    phi.prior = "fixed", phi = range_m,
    tausq.rel.prior = "fixed", tausq.rel = geor_noise[1]
  ),
  output = geoR::output.control(
    n.posterior = 0, n.predictive = 0,
    moments = TRUE, simulations.predictive = FALSE,
    signal = TRUE, messages = FALSE
  )
)

# geoR returns deterministic latent-process moments here. lambda = 1 avoids
# logging the already transformed residuals a second time.
geor_results <- stock_results(
  raw_pixels,
  raw_prior$mean + sqrt(raw_prior$variance) * as.vector(geor_fit$predictive$mean),
  raw_prior$variance * as.vector(geor_fit$predictive$variance)
)
save_analysis(geor_results, "raw_cores_geoR")

# Check both package engines directly, even without handwritten references.
geor_gstat_comparison <- do.call(rbind, lapply(
  c("posterior_log_mean", "posterior_log_var", "posterior_mean", "posterior_sd",
    "lower95", "upper95"),
  function(quantity) {
    delta <- geor_results[[quantity]] - raw_results[[quantity]]
    data.frame(
      quantity = quantity, cells = nrow(raw_pixels),
      max_absolute_difference = max(abs(delta)),
      rmse_difference = sqrt(mean(delta^2)),
      within_tolerance = all(abs(delta) <= 1e-8 + 1e-8 * abs(raw_results[[quantity]]))
    )
  }
))
write_csv(geor_gstat_comparison, file.path(output_folder, "geoR_vs_gstat.csv"))
if (!all(geor_gstat_comparison$within_tolerance)) {
  stop("geoR and gstat differ: inspect geoR_vs_gstat.csv.")
}


#==========================================
#----- 2. Area-based: mvtnorm --------------
#==========================================

# This helper delegates conditioning to mvtnorm instead of implementing
# posterior precision weights. It is used independently for each map alternative.
condition_map <- function(d, error_sd) {
  stopifnot(length(error_sd) == 1, is.finite(error_sd), error_sd > 0)
  p <- log_prior(d$prior_mean, d$prior_sd)
  error_var <- log1p((error_sd / d$field_estimate)^2)

  # Z is the standardised latent stock; U is a noisy observation of Z.
  # Joint covariance of (U, Z): [[1 + q, 1], [1, 1]], q = error_var/prior_var.
  q <- error_var / p$variance
  joint_covariance <- mvtnorm::syMatrices(
    rbind(1 + q, rep(1, length(q)), rep(1, length(q))), diag = TRUE
  )
  conditional <- mvtnorm::cond_mvnorm(
    chol = chol(joint_covariance), which_given = 1L,
    given = matrix((log(d$field_estimate) - p$mean) / sqrt(p$variance), nrow = 1L)
  )
  result <- stock_results(
    d,
    p$mean + sqrt(p$variance) * as.numeric(conditional$mean),
    p$variance * as.numeric(mvtnorm::Tcrossprod(conditional$chol, diag_only = TRUE))
  )
  result$working_error_sd <- error_sd
  result
}

area_sd <- sd(cores$stock)
area_pixels <- pixel_table(area_grid)
stopifnot(all(abs(area_pixels$field_estimate - mean(cores$stock)) < 1e-5))
area_results <- condition_map(area_pixels, area_sd)
save_analysis(area_results, "area_based")


#==========================================
#----- 3. IDW: gstat + mvtnorm -------------
#==========================================

idw_parameters <- read_csv(input_files[["idw_parameters"]], show_col_types = FALSE)
idw_parameters <- idw_parameters[idw_parameters$Variable == "stock_0_100", , drop = FALSE]
stopifnot(nrow(idw_parameters) == 1, idw_parameters$N_Cores == nrow(cores))
idw_points <- sf::st_as_sf(core_points)
heldout <- numeric(nrow(cores))
for (plot_id in unique(cores$plot_id)) {
  test <- cores$plot_id == plot_id
  heldout[test] <- gstat::idw(
    stock ~ 1, locations = idw_points[!test, ], newdata = idw_points[test, ],
    idp = idw_parameters$idp, nmax = min(idw_parameters$nmax, sum(!test)), debug.level = 0
  )$var1.pred
}
idw_checks <- data.frame(core_id = cores$core_id, plot_id = cores$plot_id,
                         observed = cores$stock, predicted = heldout, error = cores$stock - heldout)
idw_sd <- sqrt(mean(idw_checks$error^2))
idw_results <- condition_map(pixel_table(idw_grid), idw_sd)
save_analysis(idw_results, "idw")
write_csv(idw_checks, file.path(output_folder, "idw/error_checks.csv"))


#==========================================
#----- Compare with saved handwritten runs -
#==========================================

package_results <- list(raw_cores = raw_results, raw_cores_geoR = geor_results,
                        area_based = area_results, idw = idw_results)
comparison_summary <- list()
stock_plot_data <- list()
difference_plot_data <- list()

if (compare_handwritten) {
  reference_files <- c(
    raw_cores = file.path(handwritten_folder, "raw_cores/raw_core_predictions.csv"),
    raw_cores_geoR = file.path(handwritten_folder, "raw_cores/raw_core_predictions.csv"),
    area_based = file.path(handwritten_folder, "area_based/predictions.csv"),
    idw = file.path(handwritten_folder, "idw/predictions.csv")
  )
  if (!all(file.exists(reference_files))) stop("Run the handwritten analyses before comparing.")
  stopifnot(compareGeom(rast(file.path(handwritten_folder, "prepared_prior.tif")), grid_prior,
                        stopOnError = FALSE))

  # Check that the raw-core baseline used the same fixed settings.
  baseline <- readRDS(file.path(handwritten_folder, "raw_cores/raw_core_results.rds"))
  baseline_settings <- if (is.null(baseline$settings)) baseline else baseline$settings
  stopifnot(
    isTRUE(all.equal(baseline_settings$range_m, range_m)),
    isTRUE(all.equal(baseline_settings$observation_log_sd, observation_log_sd)),
    isTRUE(all.equal(baseline$cores$stock, cores$stock)),
    identical(as.character(baseline$cores$core_id), cores$core_id)
  )

  metrics <- c("prior_mean", "prior_sd", "posterior_mean", "posterior_sd",
               "lower95", "upper95", "change_from_prior", "log_variance_reduction")
  for (method in names(package_results)) {
    package <- package_results[[method]]
    written <- read_csv(reference_files[[method]], show_col_types = FALSE)
    stopifnot(!anyDuplicated(written$cell), setequal(package$cell, written$cell))
    written <- written[match(package$cell, written$cell), ]
    stopifnot(all(abs(package$x - written$x) < 1e-6), all(abs(package$y - written$y) < 1e-6))
    differences <- package[, c("cell", "x", "y")]
    for (metric in metrics) {
      delta <- package[[metric]] - written[[metric]]
      pass <- all(is.finite(delta) & abs(delta) <= 1e-8 + 1e-8 * abs(written[[metric]]))
      comparison_summary[[length(comparison_summary) + 1L]] <- data.frame(
        method = method, quantity = metric, cells = nrow(package),
        max_absolute_difference = max(abs(delta)), rmse_difference = sqrt(mean(delta^2)),
        mean_difference = mean(delta), within_tolerance = pass
      )
      differences[[metric]] <- delta
    }
    write_csv(differences, file.path(output_folder, method, "differences_from_handwritten.csv"))
    stock_plot_data[[method]] <- rbind(
      data.frame(x = package$x, y = package$y, stock = written$posterior_mean, method, engine = "Handwritten"),
      data.frame(x = package$x, y = package$y, stock = package$posterior_mean, method, engine = "Package")
    )
    delta <- package$posterior_mean - written$posterior_mean
    difference_plot_data[[method]] <- data.frame(
      x = package$x, y = package$y, method,
      difference = ifelse(abs(delta) < 1e-8, 0, delta)
    )
  }
  comparison_summary <- do.call(rbind, comparison_summary)
  write_csv(comparison_summary, file.path(output_folder, "comparison_summary.csv"))

  stock_data <- do.call(rbind, stock_plot_data)
  stock_data$method <- factor(stock_data$method, levels = c("raw_cores", "raw_cores_geoR", "area_based", "idw"),
                              labels = c("Raw cores: gstat", "Raw cores: geoR", "Area-based", "IDW"))
  comparison_plot <- ggplot(stock_data, aes(x, y, fill = stock)) +
    geom_raster() + facet_grid(method ~ engine) + coord_equal() +
    scale_fill_viridis_c(name = "kg C/m²") + theme_minimal() +
    labs(title = "Handwritten and package posterior means", x = "Easting (m)", y = "Northing (m)")
  ggsave(file.path(output_folder, "comparison_maps.png"), comparison_plot, width = 10, height = 16, dpi = 200)

  difference_data <- do.call(rbind, difference_plot_data)
  difference_data$method <- factor(difference_data$method, levels = c("raw_cores", "raw_cores_geoR", "area_based", "idw"),
                                   labels = c("Raw cores: gstat", "Raw cores: geoR", "Area-based", "IDW"))
  limit <- max(1e-8, max(abs(difference_data$difference)))
  difference_plot <- ggplot(difference_data, aes(x, y, fill = difference)) +
    geom_raster() + facet_wrap(~method, ncol = 2) + coord_equal() + theme_minimal() +
    theme(panel.background = element_rect(fill = "grey92", colour = NA)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-limit, limit),
                         name = "kg C/m²") +
    labs(title = "Package minus handwritten means", subtitle = "White cells differ by less than 1e-8; grey areas are unmapped",
         x = "Easting (m)", y = "Northing (m)")
  ggsave(file.path(output_folder, "difference_maps.png"), difference_plot, width = 10, height = 9, dpi = 200)
  print(comparison_summary)
  if (!all(comparison_summary$within_tolerance)) stop("Differences exceed tolerance; inspect comparison_summary.csv.")
}


#==========================================
#----- Summaries, provenance and references -
#==========================================

model_summary <- do.call(rbind, lapply(names(package_results), function(method) {
  d <- package_results[[method]]
  data.frame(method = method, prediction_cells = nrow(d),
             mean_prior_stock_kg_m2 = mean(d$prior_mean),
             mean_updated_stock_kg_m2 = mean(d$posterior_mean),
             minimum_stock_kg_m2 = min(d$posterior_mean), maximum_stock_kg_m2 = max(d$posterior_mean),
             mean_pointwise_sd_kg_m2 = mean(d$posterior_sd))
}))
write_csv(model_summary, file.path(output_folder, "model_summary.csv"))
write_csv(data.frame(
  method = names(package_results), package = c("gstat", "geoR", "mvtnorm", "gstat + mvtnorm"),
  error_source = c("Fixed observation log-SD", "Fixed observation log-SD",
                   "Between-core SD", "Held-out plot RMSE; IDW settings fixed"),
  working_error_value = c(observation_log_sd, observation_log_sd, area_sd, idw_sd),
  error_units = c("log-stock SD", "log-stock SD", "kg C/m²", "kg C/m²"),
  range_m = c(range_m, range_m, NA, NA),
  inference = c("Spatial Gaussian conditioning", "Spatial Gaussian conditioning; fixed parameters",
                "Approximate marginal map fusion", "Approximate marginal map fusion"),
  intervals = "Pointwise latent-stock intervals conditional on fixed assumptions"
), file.path(output_folder, "model_settings.csv"))

fingerprints <- c(input_files, if (compare_handwritten) c(
  reference_files,
  handwritten_grid = file.path(handwritten_folder, "prepared_prior.tif"),
  handwritten_raw_settings = file.path(handwritten_folder, "raw_cores/raw_core_results.rds")
) else character())
write_csv(data.frame(input = names(fingerprints), path = unname(fingerprints),
                     md5 = unname(tools::md5sum(fingerprints))), file.path(output_folder, "input_manifest.csv"))
write_csv(data.frame(package = c("terra", "gstat", "geoR", "mvtnorm", "sp", "sf", "readr", "ggplot2"),
                     version = vapply(c("terra", "gstat", "geoR", "mvtnorm", "sp", "sf", "readr", "ggplot2"),
                                      function(x) as.character(packageVersion(x)), character(1))),
          file.path(output_folder, "package_versions.csv"))
writeLines(c(capture.output(citation("gstat")), "", capture.output(citation("geoR")), "",
             capture.output(citation("mvtnorm")), "",
             "geoR Bayesian kriging: https://search.r-project.org/CRAN/refmans/geoR/html/krige.bayes.html",
             "mvtnorm package background: https://journal.r-project.org/articles/RJ-2009-001/",
             "Conditional-normal API: https://search.r-project.org/CRAN/refmans/mvtnorm/html/margcond.html"),
           file.path(output_folder, "package_references.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_folder, "session_info.txt"))

cat("\nPackage analyses completed. Outputs:", output_folder, "\n")
print(model_summary)
