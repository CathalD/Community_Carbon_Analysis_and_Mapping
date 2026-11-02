# ==========================================
# ------ Bayesian Spatial Modelling -----
# =============================================

# ===================
# Step 1 - Settings
# ===================
output_folder <- "Draft1/outputs"  # Change for each depth or sensitivity run

# Inputs from your existing workshop folder.
core_file <- read_csv(
  "Draft1/data/mapping_data.csv",
  show_col_types = FALSE
)
head(core_file)
field_points <- st_as_sf(
  core_file, 
  coords = c("Longitude", "Latitude"), 
  crs = 4326
)


aoi_file <- vect(st_read("Draft1/outputs/aoi_polygon.gpkg"))

# --------------------------------
# Spatial Output from Script 2
# ------------------------------------
library(sf)
st_layers("Draft1/outputs/area_based.gpkg")
st_layers("Draft1/outputs/aoi_polygon.gpkg")
# Step 1 - Load current maps:
aoi <- vect(st_read("Draft1/outputs/aoi_polygon.gpkg"))
# Area-based extrapolation
ab_0_30 <- st_read(
  area_based_gpkg,
  layer = "area_based_0_30"
)
ab_0_100 <- st_read(
  area_based_gpkg,
  layer = "area_based_0_100"
)
# IDW
idw_0_30 <- rast("Draft1/outputs/IDW_0_30.tif")
idw_0_100 <- rast("Draft1/outputs/IDW_0_100.tif")



# ------------------------
#    Import prior data
# -------------------------


prior_map <- rast("Draft1/outputs/prior_MM_0_100.tif")
prior_uncertainty <- rast("Draft1/outputs/prior_MM_sd.tif")
prior_var <- prior_sd^2



# Double check and and match crs for field, modelle, priors

if (st_crs(field_points) != st_crs(prior_map)) {
  field_points <- st_transform(field_points, st_crs(prior_map))
}

# ========================
# ------
#========================

# ---- FIELD ----
field_mean <- mean(
  core_file$`C Stock 0-100cm`,
  na.rm = TRUE
)

field_sd <- sd(
  core_file$`C Stock 0-100cm`,
  na.rm = TRUE
)

field_ci95 <- 1.96 * field_sd / sqrt(nrow(core_file))

field_carbon <- rast(idw_0_100)
values(field_carbon) <- field_mean

field_uncertainty <- rast(idw_0_100)
values(field_uncertainty) <- field_ci95

field_cv <- (field_uncertainty / field_carbon) * 100


# ---- AREA-BASED ----
ab_raster <- rasterize(
  vect(ab_0_100),
  idw_0_100,
  field = "mean_stock_kg_m2"
)

ab_uncertainty <- rast(ab_raster)
values(ab_uncertainty) <- field_ci95

ab_cv <- (ab_uncertainty / ab_raster) * 100


# ---- IDW ----

idw_uncertainty <- rast(idw_0_100)
values(idw_uncertainty) <- 

idw_cv <- (idw_uncertainty / idw_0_100) * 100


# ---- PRIOR ----
prior_cv <- (prior_uncertainty / prior_map) * 100


# ==========================================
# Summary Table
# ==========================================

summary_table <- data.frame(
  Method = c(
    "Field",
    "Area-Based",
    "IDW",
    "Prior"
  ),
  Mean = c(
    field_mean,
    global(ab_raster, "mean", na.rm = TRUE)[1,1],
    global(idw_0_100, "mean", na.rm = TRUE)[1,1],
    global(prior_map, "mean", na.rm = TRUE)[1,1]
  ),
  Uncertainty = c(
    field_ci95,
    field_ci95,
    rmse_100,
    global(prior_uncertainty, "mean", na.rm = TRUE)[1,1]
  ),
  CV_Percent = c(
    global(field_cv, "mean", na.rm = TRUE)[1,1],
    global(ab_cv, "mean", na.rm = TRUE)[1,1],
    global(idw_cv, "mean", na.rm = TRUE)[1,1],
    global(prior_cv, "mean", na.rm = TRUE)[1,1]
  )
)

summary_table














target_depth_cm <- 100                # 30, 50 or 100; change both prior files too
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
depth_label <- paste0("0-", target_depth_cm, " cm below ground surface")
prior_source <- "Sothe et al (2022)"



# Configuration:
uncertainty_type <- "sd"              # sd, variance, pi90_width, pi95_width
mean_to_kg_m2 <- 1                    # e.g. 0.1 for t C/ha -> kg C/m2
uncertainty_to_kg_m2 <- 1              # Same linear unit conversion; see README
prior_sd_multiplier <- 1              # Sensitivity factor after conversion
map_crs <- 3978                       # Canada Atlas Lambert, metres
cell_size_m <- 1000                   # Prediction spacing, not proven accuracy
range_m <- 5000                       # Exponential e-folding distance
observation_log_sd <- 0.15             # Assumed log-scale observation/support SD
chunk_size <- 5000L

# Install once if you need to
library(sf)
library(terra)
library(ggplot2)
library(patchwork)

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# Step 3 - Read the same core data as Scripts 1 and 2
# =======================================================
raw_cores <- (core_file)

aoi <- (aoi_file)

stopifnot(all(c("Core ID", "Plot ID", "Longitude", "Latitude", stock_column) %in% names(raw_cores)))
cores <- data.frame(core_id = raw_cores$`Core ID`, group_id = raw_cores$`Plot ID`,
                    longitude = raw_cores$Longitude, latitude = raw_cores$Latitude,
                    stock = raw_cores[[stock_column]])
cores$source_notes <- if ("Notes" %in% names(raw_cores)) raw_cores$Notes else ""
print(cores[, c("core_id", "stock")])



# ===================
# Step 5 - Read the prior maps and convert uncertainty to stock SD
# ===================

# Read only the AOI window when the input maps cover all of Canada.
prior_extent <- ext(project((aoi), crs(prior_mean_file)))

prior_mean <- crop(prior_mean_file, prior_extent, snap = "out")

prior_uncertainty <- crop(prior_uncertainty_file, prior_extent, snap = "out")

cat("Prior:", prior_source, "\nTarget:", depth_label, "\n")
prior_mean <- prior_mean * mean_to_kg_m2
prior_sd <- switch(uncertainty_type,
                   sd = prior_uncertainty * uncertainty_to_kg_m2,
                   variance = sqrt(prior_uncertainty) * uncertainty_to_kg_m2,
                   pi90_width = prior_uncertainty * uncertainty_to_kg_m2 / (2 * qnorm(0.95)),
                   pi95_width = prior_uncertainty * uncertainty_to_kg_m2 / (2 * qnorm(0.975)),
                   stop("Unknown uncertainty_type. Relative indices cannot be used as stock SD."))
prior_sd <- prior_sd * prior_sd_multiplier



# Sample the native prior at the cores, before making a coarser display grid.
native_prior <- c(prior_mean, prior_sd)
names(native_prior) <- c("prior_mean", "prior_sd")
at_cores <- extract(native_prior, vect(core_sf), ID = FALSE)
valid_prior <- is.finite(at_cores$prior_mean) & at_cores$prior_mean > 0 &
  is.finite(at_cores$prior_sd) & at_cores$prior_sd > 0
missing_prior <- cores[!valid_prior, ]
missing_prior$reason <- rep("prior mean or SD missing/nonpositive", nrow(missing_prior))
excluded <- rbind(excluded, missing_prior)
write.csv(excluded, file.path(output_folder, "excluded_cores.csv"), row.names = FALSE)
cores <- cores[valid_prior, ]; core_sf <- core_sf[valid_prior, ]
at_cores <- at_cores[valid_prior, , drop = FALSE]
if (length(unique(cores$group_id)) < 2) stop("Need at least two sampling groups for the held-out check.")
if (nrow(cores) < 5) warning("Few cores: validation is descriptive, not evidence of reliable accuracy.")
cores <- cbind(cores, at_cores)
write.csv(cores, file.path(output_folder, "included_cores.csv"), row.names = FALSE)
print(cores[, c("core_id", "stock", "prior_mean", "prior_sd")])

# ===================
# Step 6 - Build a prediction grid
# ===================
box <- st_bbox(aoi)
template <- rast(xmin = floor(box["xmin"] / cell_size_m) * cell_size_m,
                 xmax = ceiling(box["xmax"] / cell_size_m) * cell_size_m,
                 ymin = floor(box["ymin"] / cell_size_m) * cell_size_m,
                 ymax = ceiling(box["ymax"] / cell_size_m) * cell_size_m,
                 resolution = cell_size_m, crs = st_crs(aoi)$wkt)
if (ncell(template) > 250000) stop("Use a coarser cell_size_m for this workshop (maximum 250,000 cells).")
# Nearest-cell sampling preserves each selected mean/SD pair. It is not block averaging.
grid_prior <- project(native_prior, template, method = "near")
grid_prior <- mask(grid_prior, vect(aoi), touches = FALSE)
grid_values <- values(grid_prior)
eligible_cells <- which(is.finite(grid_values[, 1]) & grid_values[, 1] > 0 &
                          is.finite(grid_values[, 2]) & grid_values[, 2] > 0)
if (!length(eligible_cells)) stop("No valid prior cells inside the AOI.")
grid_xy <- xyFromCell(template, eligible_cells)
core_xy <- st_coordinates(core_sf)
print(ggplot() + geom_sf(data = aoi, fill = "grey95") +
        geom_sf(data = core_sf, aes(colour = stock), size = 3) +
        scale_colour_viridis_c(name = "kg C/m²") + theme_minimal())

# ===================
# Step 7 - Express the prior distribution on the log scale
# ===================
# These equations preserve the original arithmetic mean and SD under a lognormal prior.
core_log_var <- log1p((cores$prior_sd / cores$prior_mean)^2)
core_log_mean <- log(cores$prior_mean) - core_log_var / 2
grid_log_var <- log1p((grid_values[eligible_cells, 2] / grid_values[eligible_cells, 1])^2)
grid_log_mean <- log(grid_values[eligible_cells, 1]) - grid_log_var / 2
observed_log <- log(cores$stock)
# log(observed stock) = log(true stock) + Normal(0, observation_log_sd^2).
# Errors are unbiased on the log scale, not the arithmetic stock scale.
observation_var <- rep(observation_log_sd^2, nrow(cores))
core_distances <- as.matrix(dist(core_xy))

# ===================
# Step 8 - Condition a Gaussian spatial prior on the observations
# ===================
# This small function is the statistical model. Cholesky solves avoid matrix inversion.
# Exponential correlation exp(-distance/range); at range the correlation is about 0.37.
condition_spatial <- function(train_mean, train_var, observed, observation_var,
                              train_distance, prediction_mean, prediction_var,
                              prediction_distance, range_m) {
  K <- outer(sqrt(train_var), sqrt(train_var)) * exp(-train_distance / range_m)
  V <- K + diag(observation_var, nrow(K))
  R <- chol(V)
  alpha <- backsolve(R, forwardsolve(t(R), observed - train_mean))
  K_new <- outer(sqrt(prediction_var), sqrt(train_var)) * exp(-prediction_distance / range_m)
  updated_mean <- prediction_mean + as.vector(K_new %*% alpha)
  projected <- forwardsolve(t(R), t(K_new))
  updated_var <- prediction_var - colSums(projected^2)
  if (any(updated_var < -1e-8 * pmax(1, prediction_var))) stop("Invalid posterior variance.")
  list(mean = updated_mean, variance = pmax(0, updated_var))
}

# ===================
# Step 9 - Predict in small batches, then return to kg C/m²
# ===================
posterior_log_mean <- posterior_log_var <- nearest_core_m <- numeric(length(eligible_cells))
for (start in seq(1, length(eligible_cells), by = chunk_size)) {
  rows <- start:min(start + chunk_size - 1, length(eligible_cells))
  distances <- sqrt(outer(grid_xy[rows, 1], core_xy[, 1], "-")^2 +
                      outer(grid_xy[rows, 2], core_xy[, 2], "-")^2)
  updated <- condition_spatial(core_log_mean, core_log_var, observed_log,
                               observation_var, core_distances, grid_log_mean[rows], grid_log_var[rows], distances, range_m)
  posterior_log_mean[rows] <- updated$mean
  posterior_log_var[rows] <- updated$variance
  nearest_core_m[rows] <- apply(distances, 1, min)
}
posterior_mean <- exp(posterior_log_mean + posterior_log_var / 2)
posterior_median <- exp(posterior_log_mean)
posterior_sd <- sqrt(expm1(posterior_log_var) * exp(2 * posterior_log_mean + posterior_log_var))
lower95 <- exp(posterior_log_mean - qnorm(0.975) * sqrt(posterior_log_var))
upper95 <- exp(posterior_log_mean + qnorm(0.975) * sqrt(posterior_log_var))

# ===================
# Step 10 - Leave one sampling group out
# ===================
# Withhold all cores in a plot together. The prior must also be independent of these cores.
validation <- list()
for (group in unique(cores$group_id)) {
  test <- which(cores$group_id == group)
  train <- which(cores$group_id != group)
  heldout <- condition_spatial(core_log_mean[train], core_log_var[train], observed_log[train],
                               observation_var[train], core_distances[train, train, drop = FALSE],
                               core_log_mean[test], core_log_var[test], core_distances[test, train, drop = FALSE], range_m)
  # A new observed core also has measurement/support variation.
  predictive_var <- heldout$variance + observation_var[test]
  validation[[length(validation) + 1]] <- data.frame(
    core_id = cores$core_id[test], group_id = group, observed = cores$stock[test],
    prior = exp(core_log_mean[test] + (core_log_var[test] + observation_var[test]) / 2),
    bayes = exp(heldout$mean + predictive_var / 2),
    lower95 = exp(heldout$mean - qnorm(0.975) * sqrt(predictive_var)),
    upper95 = exp(heldout$mean + qnorm(0.975) * sqrt(predictive_var)))
}
validation <- do.call(rbind, validation)
write.csv(validation, file.path(output_folder, "heldout_predictions.csv"), row.names = FALSE)
scores <- do.call(rbind, lapply(c("prior", "bayes"), function(method) {
  error <- validation[[method]] - validation$observed
  data.frame(method = method, n = nrow(validation), MAE = mean(abs(error)),
             RMSE = sqrt(mean(error^2)), bias = mean(error))
}))
print(scores)
cat("Held-out 95% predictive coverage:",
    mean(validation$observed >= validation$lower95 & validation$observed <= validation$upper95), "\n")
write.csv(scores, file.path(output_folder, "validation_scores.csv"), row.names = FALSE)
# Repeat with 2.5, 5 and 10 km ranges in separate output folders to assess sensitivity.

# ===================
# Step 11 - Save the maps
# ===================
result <- cbind(prior_mean = grid_values[eligible_cells, 1],
                prior_sd = grid_values[eligible_cells, 2], posterior_mean = posterior_mean,
                posterior_median = posterior_median, posterior_sd = posterior_sd,
                lower95 = lower95, upper95 = upper95, interval_width95 = upper95 - lower95,
                change_from_prior = posterior_mean - grid_values[eligible_cells, 1],
                log_variance_reduction_fraction = 1 - posterior_log_var / grid_log_var,
                nearest_core_m = nearest_core_m)
if (any(!is.finite(result))) stop("Non-finite output; inspect input distributions and assumptions.")
maps <- rast(template, nlyrs = ncol(result))
names(maps) <- colnames(result)
all_values <- matrix(NA_real_, ncell(template), ncol(result))
all_values[eligible_cells, ] <- result
values(maps) <- all_values
writeRaster(maps, file.path(output_folder, "bayesian_carbon_maps.tif"), overwrite = TRUE)

# ===================
# Step 12 - Make a workshop figure
# ===================
plot_data <- data.frame(grid_xy, result)
names(plot_data)[1:2] <- c("x", "y")
stock_limits <- range(c(result[, "prior_mean"], posterior_mean))
make_map <- function(variable, title, limits = NULL, diverging = FALSE) {
  p <- ggplot(plot_data, aes(x, y, fill = .data[[variable]])) +
    geom_tile(width = cell_size_m, height = cell_size_m) +
    geom_sf(data = aoi, inherit.aes = FALSE, fill = NA, linewidth = 0.3) +
    geom_sf(data = core_sf, inherit.aes = FALSE, shape = 21, fill = "white", size = 2) +
    coord_sf(crs = st_crs(aoi)) + labs(title = title, fill = "kg C/m²", x = NULL, y = NULL) +
    theme_minimal(base_size = 10)
  if (diverging) p + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)
  else p + scale_fill_viridis_c(limits = limits)
}
panel <- (make_map("prior_mean", "Prior mean", stock_limits) |
            make_map("posterior_mean", "Updated mean", stock_limits)) /
  (make_map("change_from_prior", "Change from prior", diverging = TRUE) |
     make_map("interval_width95", "Updated 95% interval width")) +
  plot_annotation(title = paste("Carbon map update:", depth_label),
                  subtitle = "Conditional spatial model; white points are included cores",
                  caption = "Input provenance is retained in included_cores.csv; see the workshop README.")
print(panel)
ggsave(file.path(output_folder, "bayesian_workshop_maps.png"), panel, width = 11, height = 9, dpi = 200)

# ===================
# Step 13 - Record assumptions, inputs and software
# ===================
settings <- list(target_depth_cm = target_depth_cm, depth_label = depth_label, prior_source = prior_source,
                 uncertainty_type = uncertainty_type, mean_to_kg_m2 = mean_to_kg_m2,
                 uncertainty_to_kg_m2 = uncertainty_to_kg_m2, prior_sd_multiplier = prior_sd_multiplier,
                 range_m = range_m,
                 observation_log_sd = observation_log_sd, map_crs = map_crs, cell_size_m = cell_size_m,
                 stock_column = stock_column,
                 depth_assumption = "Use Script 1 intervals ending at or before target; retain shallower full profiles",
                 inference = "Exact conditional log-Gaussian update; covariance assumptions fixed",
                 spatial_support = "Point predictions at grid centres; nearest native prior cell",
                 observation_model = "Log-unbiased multiplicative error; predictive means include observation variance")
saveRDS(list(settings = settings, cores = cores, log_mean = posterior_log_mean,
             log_variance = posterior_log_var, eligible_cells = eligible_cells,
             validation = validation, scores = scores), file.path(output_folder, "workshop_results.rds"))
capture.output(dput(settings), file = file.path(output_folder, "settings.txt"))
write.csv(data.frame(path = input_paths,
                     md5 = unname(tools::md5sum(input_paths))), file.path(output_folder, "input_checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_folder, "session_info.txt"))
capture.output(for (package in c("base", "sf", "terra", "ggplot2", "patchwork")) {
  cat("\nPACKAGE:", package, "\n"); print(citation(package))
}, file = file.path(output_folder, "package_citations.txt"))
cat("\nFinished. Outputs:", normalizePath(output_folder), "\n")
cat("Intervals are conditional on the recorded assumptions. See README before interpretation.\n")