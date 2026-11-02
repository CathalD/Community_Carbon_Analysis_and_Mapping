#==========================================
#----- Script 5: Bayesian carbon mapping ---
#==========================================

library(terra)
library(readr)
library(dplyr)
library(ggplot2)
terra::projNetwork(FALSE)


#==========================================
#----- 5A: Settings and inputs -------------
#==========================================

target_depth_cm <- 100

stock_column <- paste0(
  "C Stock 0-", target_depth_cm, "cm"
)

output_folder <- file.path(
  "outputs/script5",
  paste0("stock_0_", target_depth_cm)
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


#----- Input files -----

input_files <- c(
  field = "outputs/script1/mapping_data.csv",
  boundary = "data/boundary/aoi.geojson",
  prior_mean = "outputs/script4/mcmaster_0_100cm.tif",
  prior_uncertainty =
    "outputs/script4/mcmaster_uncertainty_0_100cm.tif",
  grid = "outputs/script2/idw_0_100cm.tif"
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}


#----- Field observations -----

raw_cores <- read_csv(
  input_files[["field"]],
  show_col_types = FALSE
)

required_columns <- c(
  "Core ID", "Plot ID",
  "Longitude", "Latitude",
  stock_column
)

stopifnot(
  all(required_columns %in% names(raw_cores))
)

cores <- data.frame(
  core_id = as.character(raw_cores[["Core ID"]]),
  plot_id = as.character(raw_cores[["Plot ID"]]),
  longitude = raw_cores[["Longitude"]],
  latitude = raw_cores[["Latitude"]],
  stock = raw_cores[[stock_column]]
)


#----- Prior maps -----

prior_mean <- rast(
  input_files[["prior_mean"]]
)

prior_uncertainty <- rast(
  input_files[["prior_uncertainty"]]
)

names(prior_mean) <- "prior_mean"
names(prior_uncertainty) <- "prior_interval_width_90"

# Prepare the conversion in 5B.
# These are interval widths, in kg C/m².
prior_uncertainty_type <- "pi90_width"


#----- Output grid and boundary -----

# Use the grid geometry established in Script 2.
template <- rast(input_files[["grid"]])

aoi <- vect(input_files[["boundary"]])
aoi <- project(aoi, crs(template))

core_points <- vect(
  cores,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

core_points <- project(
  core_points,
  crs(template)
)


#----- Check what was loaded -----

print(cores[, c("core_id", "plot_id", "stock")])

cat("\nTarget depth: 0–", target_depth_cm, " cm\n", sep = "")
cat("Field cores:", nrow(cores), "\n")
cat("Sampling plots:", length(unique(cores$plot_id)), "\n")
cat("Grid resolution:", res(template), "metres\n")
cat("Stock units: kg C/m²\n")



#==========================================
#----- 5B: Prepare cores and prior ---------
#==========================================

#----- Convert interval width to working SD -----

stopifnot(prior_uncertainty_type == "pi90_width")

if (!compareGeom(
  prior_mean,
  prior_uncertainty,
  stopOnError = FALSE
)) {
  stop("The prior mean and uncertainty maps must share a grid.")
}

# Approximate conversion assuming a normal-shaped interval.
# This working SD is not a calibrated measurement-error estimate.

prior_sd <- prior_uncertainty / (2 * qnorm(0.95))

native_prior <- c(prior_mean, prior_sd)

names(native_prior) <- c(
  "prior_mean",
  "prior_sd"
)


#----- Extract prior values at each core -----

at_cores <- terra::extract(
  native_prior,
  project(core_points, crs(native_prior)),
  method = "simple",
  ID = FALSE
)

inside_aoi <- rowSums(
  relate(core_points, aoi, "intersects")
) > 0


#----- Check each observation -----

# Positive stocks are required for the log-scale model.
# Duplicate core IDs are flagged to avoid counting cores twice.

core_checks <- bind_cols(cores, at_cores) %>%
  mutate(
    core_id = trimws(core_id),
    plot_id = trimws(plot_id),
    
    valid_stock = is.finite(stock) & stock > 0,
    
    valid_ids =
      !is.na(core_id) & nzchar(core_id) &
      !is.na(plot_id) & nzchar(plot_id),
    
    unique_core_id = !(
      duplicated(core_id) |
        duplicated(core_id, fromLast = TRUE)
    ),
    
    within_aoi = inside_aoi,
    
    valid_prior =
      is.finite(prior_mean) & prior_mean > 0 &
      is.finite(prior_sd) & prior_sd > 0,
    
    included =
      valid_stock &
      valid_ids &
      unique_core_id &
      within_aoi &
      valid_prior
  )

write_csv(
  core_checks,
  file.path(output_folder, "core_checks.csv")
)

write_csv(
  filter(core_checks, !included),
  file.path(output_folder, "excluded_cores.csv")
)


#----- Retain usable cores -----

cores <- core_checks %>%
  filter(included) %>%
  select(
    core_id,
    plot_id,
    longitude,
    latitude,
    stock,
    prior_mean,
    prior_sd
  )

core_points <- core_points[core_checks$included, ]

write_csv(
  cores,
  file.path(output_folder, "included_cores.csv")
)

if (n_distinct(cores$plot_id) < 2) {
  stop(
    "At least two usable sampling plots are needed ",
    "for the held-out evaluation."
  )
}


#----- Prepare prior values on the output grid -----

grid_prior <- native_prior

if (!compareGeom(
  native_prior[[1]],
  template,
  stopOnError = FALSE
)) {
  
  # Sample the mean and SD together using nearest cells.
  # This preserves paired values; it does not estimate
  # uncertainty for an aggregated 1 km block.
  
  grid_prior <- project(
    native_prior,
    template,
    method = "near"
  )
}

grid_prior <- mask(
  grid_prior,
  project(aoi, crs(grid_prior)),
  touches = FALSE
)

valid_grid <- (
  is.finite(grid_prior[["prior_mean"]]) &
    grid_prior[["prior_mean"]] > 0 &
    is.finite(grid_prior[["prior_sd"]]) &
    grid_prior[["prior_sd"]] > 0
)

grid_prior <- mask(
  grid_prior,
  valid_grid,
  maskvalues = 0
)

if (
  global(valid_grid, "sum", na.rm = TRUE)[1, 1] == 0
) {
  stop("No usable prior cells remain inside the boundary.")
}

writeRaster(
  grid_prior,
  file.path(output_folder, "prepared_prior.tif"),
  overwrite = TRUE
)


#----- Record preparation assumptions -----

writeLines(
  c(
    "Stock units: kg C/m2.",
    "Prior uncertainty input: Q95 minus Q05.",
    "Working SD: interval width / (2 * qnorm(0.95)).",
    "The SD conversion assumes a normal-shaped interval.",
    "The supplied prior rasters have already been resampled.",
    "Working prior spread is approximate for this demonstration.",
    "Grid values do not establish uncertainty for block averages.",
    "Field stocks retain the definitions established in Script 1."
  ),
  file.path(output_folder, "preparation_notes.txt")
)


#----- Review the prepared observations -----

print(
  cores[, c(
    "core_id",
    "stock",
    "prior_mean",
    "prior_sd"
  )]
)

cat("\nIncluded cores:", nrow(cores), "\n")
cat("Excluded cores:", sum(!core_checks$included), "\n")
cat("Sampling plots:", n_distinct(cores$plot_id), "\n")


#==========================================
#----- 5C: Spatial Bayesian update ---------
#==========================================

#----- Working model settings -----

range_m <- 2000
observation_log_sd <- 0.15
chunk_size <- 5000L

stopifnot(
  range_m > 0,
  observation_log_sd > 0,
  nrow(cores) == nrow(core_points)
)

# Coordinates use the workshop's projected grid in metres.
if (is.lonlat(grid_prior)) {
  stop("The prediction grid must use projected coordinates in metres.")
}


#----- Prepare prediction locations -----

grid_values <- values(grid_prior, mat = TRUE)

eligible_cells <- which(
  is.finite(grid_values[, 1]) &
    grid_values[, 1] > 0 &
    is.finite(grid_values[, 2]) &
    grid_values[, 2] > 0
)

grid_xy <- xyFromCell(
  grid_prior,
  eligible_cells
)

core_xy <- crds(
  project(core_points, crs(grid_prior))
)


#----- Express the prior on the log scale -----

# These equations preserve the working arithmetic mean
# and SD when representing the prior as lognormal.

core_log_var <- log1p(
  (cores$prior_sd / cores$prior_mean)^2
)

core_log_mean <- log(cores$prior_mean) - core_log_var / 2

grid_log_var <- log1p(
  (
    grid_values[eligible_cells, 2] /
      grid_values[eligible_cells, 1]
  )^2
)

grid_log_mean <- (
  log(grid_values[eligible_cells, 1]) -
    grid_log_var / 2
)

observed_log <- log(cores$stock)

# Observation model:
# log(observed stock) = log(underlying stock) + error
# error follows Normal(0, observation_log_sd^2).
# This assumes unbiased errors on the log scale.

observation_var <- rep(
  observation_log_sd^2,
  nrow(cores)
)


#----- Calculate distances between locations -----

distance_matrix <- function(a, b) {
  
  sqrt(
    outer(a[, 1], b[, 1], "-")^2 +
      outer(a[, 2], b[, 2], "-")^2
  )
}


#----- Define the Bayesian update -----

# Exponential spatial correlation: exp(-distance / range).
# This function will also be used for held-out evaluation.

spatial_update <- function(
    train_xy,
    train_mean,
    train_var,
    observed,
    observation_var,
    prediction_xy,
    prediction_mean,
    prediction_var,
    range_m) {
  
  # Prior covariance between training cores
  K <- (
    outer(sqrt(train_var), sqrt(train_var)) *
      exp(-distance_matrix(train_xy, train_xy) / range_m)
  )
  
  # Add observation/support variation
  V <- K + diag(
    observation_var,
    nrow = length(train_mean)
  )
  
  # Solve the model without explicitly inverting V
  R <- chol(V)
  
  alpha <- backsolve(
    R,
    forwardsolve(
      t(R),
      observed - train_mean
    )
  )
  
  # Covariance between prediction locations and cores
  K_new <- (
    outer(sqrt(prediction_var), sqrt(train_var)) *
      exp(
        -distance_matrix(prediction_xy, train_xy) / range_m
      )
  )
  
  updated_mean <- (
    prediction_mean +
      as.vector(K_new %*% alpha)
  )
  
  projected <- forwardsolve(
    t(R),
    t(K_new)
  )
  
  updated_var <- (
    prediction_var -
      colSums(projected^2)
  )
  
  # Allow only negligible floating-point differences below zero
  if (any(
    updated_var < -1e-8 * pmax(1, prediction_var)
  )) {
    stop("Invalid posterior variance; inspect the model inputs.")
  }
  
  list(
    mean = updated_mean,
    variance = pmax(0, updated_var)
  )
}


#----- Update the prediction grid in batches -----

posterior_log_mean <- numeric(length(eligible_cells))
posterior_log_var <- numeric(length(eligible_cells))

for (start in seq(
  1,
  length(eligible_cells),
  by = chunk_size
)) {
  
  rows <- start:min(
    start + chunk_size - 1,
    length(eligible_cells)
  )
  
  updated <- spatial_update(
    train_xy = core_xy,
    train_mean = core_log_mean,
    train_var = core_log_var,
    observed = observed_log,
    observation_var = observation_var,
    prediction_xy = grid_xy[rows, , drop = FALSE],
    prediction_mean = grid_log_mean[rows],
    prediction_var = grid_log_var[rows],
    range_m = range_m
  )
  
  posterior_log_mean[rows] <- updated$mean
  posterior_log_var[rows] <- updated$variance
}


#----- Return results to kg C/m² -----

posterior_mean <- exp(
  posterior_log_mean + posterior_log_var / 2
)

posterior_median <- exp(posterior_log_mean)

posterior_sd <- sqrt(
  expm1(posterior_log_var) *
    exp(2 * posterior_log_mean + posterior_log_var)
)

# Pointwise 95% credible intervals for underlying stock.
# Future-observation intervals will also include
# observation/support variation in section 5E.

lower95 <- exp(
  posterior_log_mean -
    qnorm(0.975) * sqrt(posterior_log_var)
)

upper95 <- exp(
  posterior_log_mean +
    qnorm(0.975) * sqrt(posterior_log_var)
)


#----- Assemble results for mapping in 5D -----

posterior_values <- data.frame(
  cell = eligible_cells,
  prior_mean = grid_values[eligible_cells, 1],
  prior_sd = grid_values[eligible_cells, 2],
  posterior_mean = posterior_mean,
  posterior_median = posterior_median,
  posterior_sd = posterior_sd,
  lower95 = lower95,
  upper95 = upper95,
  interval_width95 = upper95 - lower95,
  change_from_prior = (
    posterior_mean -
      grid_values[eligible_cells, 1]
  ),
  log_variance_reduction = (
    1 - posterior_log_var / grid_log_var
  )
)

if (!all(is.finite(as.matrix(posterior_values)))) {
  stop("Non-finite predictions; inspect the model inputs.")
}


#----- Retain settings for the final export -----

model_settings <- list(
  target_depth_cm = target_depth_cm,
  range_m = range_m,
  observation_log_sd = observation_log_sd,
  prior_uncertainty_type = prior_uncertainty_type,
  model = "Conditional log-Gaussian spatial update",
  spatial_correlation = "Exponential",
  spatial_support = "Point predictions at grid centres",
  interval_type = "Pointwise latent-stock credible intervals",
  parameter_treatment = "Fixed assumptions; not estimated",
  prior_sd_conversion = "90% interval width / (2 * qnorm(0.95))"
)


#----- Check the results -----

cat("\nCores used:", nrow(cores), "\n")
cat("Prediction cells:", nrow(posterior_values), "\n")

cat(
  "Updated mean stock range:",
  round(range(posterior_mean), 2),
  "kg C/m²\n"
)

cat(
  "Updated stock SD range:",
  round(range(posterior_sd), 2),
  "kg C/m²\n"
)

#==========================================
#----- 5D: Export and map the update -------
#==========================================

#----- Create the output raster layers -----

map_columns <- setdiff(names(posterior_values), "cell")

bayesian_maps <- rast(
  grid_prior[[1]],
  nlyrs = length(map_columns)
)

names(bayesian_maps) <- map_columns

map_values <- matrix(
  NA_real_,
  nrow = ncell(bayesian_maps),
  ncol = length(map_columns)
)

map_values[posterior_values$cell, ] <- as.matrix(
  posterior_values[, map_columns]
)

values(bayesian_maps) <- map_values

writeRaster(
  bayesian_maps,
  file.path(output_folder, "bayesian_carbon_maps.tif"),
  overwrite = TRUE
)

write_csv(
  posterior_values,
  file.path(output_folder, "grid_predictions.csv")
)

write_csv(
  data.frame(
    setting = names(model_settings),
    value = as.character(unlist(model_settings))
  ),
  file.path(output_folder, "model_settings.csv")
)


#----- Common colour scales -----

stock_limits <- range(
  c(
    posterior_values$prior_mean,
    posterior_values$posterior_mean
  )
)

sd_limits <- range(
  c(
    posterior_values$prior_sd,
    posterior_values$posterior_sd
  )
)

# Keep plotting limits usable if every value is identical
if (diff(stock_limits) == 0) {
  stock_limits <- stock_limits + c(-0.5, 0.5)
}

if (diff(sd_limits) == 0) {
  sd_limits <- sd_limits + c(-0.5, 0.5)
}

difference_limit <- max(
  abs(posterior_values$change_from_prior)
)

if (difference_limit == 0) difference_limit <- 1

stock_colours <- hcl.colors(100, "YlGnBu")

difference_colours <- colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(101)

map_boundary <- project(aoi, crs(bayesian_maps))
map_points <- project(core_points, crs(bayesian_maps))


#----- Small plotting helper -----

plot_bayes_layer <- function(
    layer, title, limits = NULL,
    colours = stock_colours) {
  
  plot(
    bayesian_maps[[layer]],
    main = title,
    range = limits,
    col = colours,
    plg = list(title = "kg C/m²"),
    mar = c(2, 2, 3, 5)
  )
  
  plot(
    map_boundary,
    add = TRUE,
    border = "grey30",
    lwd = 0.7
  )
  
  points(
    map_points,
    pch = 21,
    bg = "black",
    col = "white",
    cex = 0.8
  )
}


#----- Save the comparison figure -----

save_bayesian_figure <- function() {
  
  png(
    file.path(output_folder, "bayesian_maps.png"),
    width = 2400,
    height = 1600,
    res = 200
  )
  
  export_device <- dev.cur()
  on.exit(dev.off(export_device))
  
  par(
    mfrow = c(2, 3),
    oma = c(0, 0, 3, 0)
  )
  
  plot_bayes_layer(
    "prior_mean", "Prior mean",
    stock_limits
  )
  
  plot_bayes_layer(
    "posterior_mean", "Updated mean",
    stock_limits
  )
  
  plot_bayes_layer(
    "change_from_prior", "Updated minus prior",
    c(-difference_limit, difference_limit),
    difference_colours
  )
  
  plot_bayes_layer(
    "prior_sd", "Prior SD (approximate)",
    sd_limits
  )
  
  plot_bayes_layer(
    "posterior_sd", "Updated SD (conditional)",
    sd_limits
  )
  
  plot_bayes_layer(
    "interval_width95", "Updated 95% interval width"
  )
  
  mtext(
    paste0(
      "Carbon stock: 0–", target_depth_cm,
      " cm | Spatial range: ", range_m / 1000,
      " km | Points: included field cores"
    ),
    outer = TRUE,
    cex = 1
  )
}

# Close unfinished plotting devices before rebuilding the figure.
graphics.off()

save_bayesian_figure()



#==========================================
#----- 5E: Held-out evaluation -------------
#==========================================

# Each sampling plot is withheld completely.
# For new field data, this evaluation assumes the published
# prior did not use those same observations.

range_values <- sort(
  unique(c(1000, range_m, 5000))
)

validation_records <- list()


#----- Predict withheld cores for each range -----

for (test_range in range_values) {
  
  for (heldout_plot in unique(cores$plot_id)) {
    
    test <- which(cores$plot_id == heldout_plot)
    train <- which(cores$plot_id != heldout_plot)
    
    heldout <- spatial_update(
      train_xy = core_xy[train, , drop = FALSE],
      train_mean = core_log_mean[train],
      train_var = core_log_var[train],
      observed = observed_log[train],
      observation_var = observation_var[train],
      
      prediction_xy = core_xy[test, , drop = FALSE],
      prediction_mean = core_log_mean[test],
      prediction_var = core_log_var[test],
      
      range_m = test_range
    )
    
    # Predictive variance for a measured core includes
    # observation/support variation.
    
    prior_predictive_var <- (
      core_log_var[test] +
        observation_var[test]
    )
    
    bayes_predictive_var <- (
      heldout$variance +
        observation_var[test]
    )
    
    # Both methods predict observations under the same
    # log-unbiased observation model.
    
    validation_records[[length(validation_records) + 1]] <-
      data.frame(
        range_m = test_range,
        core_id = cores$core_id[test],
        plot_id = heldout_plot,
        observed = cores$stock[test],
        
        prior_prediction = exp(
          core_log_mean[test] + prior_predictive_var / 2
        ),
        
        prior_lower95 = exp(
          core_log_mean[test] -
            qnorm(0.975) * sqrt(prior_predictive_var)
        ),
        
        prior_upper95 = exp(
          core_log_mean[test] +
            qnorm(0.975) * sqrt(prior_predictive_var)
        ),
        
        bayes_prediction = exp(
          heldout$mean + bayes_predictive_var / 2
        ),
        
        bayes_lower95 = exp(
          heldout$mean -
            qnorm(0.975) * sqrt(bayes_predictive_var)
        ),
        
        bayes_upper95 = exp(
          heldout$mean +
            qnorm(0.975) * sqrt(bayes_predictive_var)
        )
      )
  }
}

validation_predictions <- bind_rows(validation_records)

write_csv(
  validation_predictions,
  file.path(output_folder, "heldout_predictions.csv")
)


#----- Put both methods into the same table -----

validation_long <- bind_rows(
  
  validation_predictions %>%
    transmute(
      range_m, core_id, plot_id, observed,
      method = "Prior",
      predicted = prior_prediction,
      lower95 = prior_lower95,
      upper95 = prior_upper95
    ),
  
  validation_predictions %>%
    transmute(
      range_m, core_id, plot_id, observed,
      method = "Bayesian",
      predicted = bayes_prediction,
      lower95 = bayes_lower95,
      upper95 = bayes_upper95
    )
)


#----- Calculate evaluation statistics -----

# Stock errors and interval widths are in kg C/m².
# Positive bias means overprediction.

validation_scores <- validation_long %>%
  group_by(range_m, method) %>%
  summarise(
    n_cores = n(),
    
    RMSE = sqrt(
      mean((predicted - observed)^2)
    ),
    
    MAE = mean(abs(predicted - observed)),
    
    bias = mean(predicted - observed),
    
    coverage95_percent = 100 * mean(
      observed >= lower95 &
        observed <= upper95
    ),
    
    mean_interval_width95 = mean(
      upper95 - lower95
    ),
    
    .groups = "drop"
  )

write_csv(
  validation_scores,
  file.path(output_folder, "validation_scores.csv")
)


#----- Report the current range separately -----

current_scores <- validation_scores[
  validation_scores$range_m == range_m,
]

current_predictions <- validation_long[
  validation_long$range_m == range_m,
]

print(current_scores)

write_csv(
  current_scores,
  file.path(output_folder, "current_validation_scores.csv")
)


#----- Observed versus held-out predictions -----

prediction_limits <- range(
  c(
    current_predictions$observed,
    current_predictions$lower95,
    current_predictions$upper95
  ),
  finite = TRUE
)

validation_plot <- ggplot(
  current_predictions,
  aes(x = observed, y = predicted)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_errorbar(
    aes(ymin = lower95, ymax = upper95),
    width = 0,
    colour = "#1F78B4",
    alpha = 0.6
  ) +
  geom_point(
    size = 2.5,
    colour = "#1F78B4"
  ) +
  facet_wrap(~ method) +
  coord_equal(
    xlim = prediction_limits,
    ylim = prediction_limits
  ) +
  labs(
    title = paste(
      "Held-out predictions: spatial range",
      range_m / 1000, "km"
    ),
    x = "Observed stock (kg C/m²)",
    y = "Predicted stock (kg C/m²)",
    caption = paste(
      "Whole sampling plots withheld.",
      "Bars: 95% predictive intervals."
    )
  ) +
  theme_minimal(base_size = 12)

print(validation_plot)

ggsave(
  file.path(output_folder, "heldout_predictions.png"),
  validation_plot,
  width = 9,
  height = 5,
  dpi = 300
)


#----- Compare alternative spatial ranges -----

score_plot_data <- bind_rows(
  
  validation_scores %>%
    transmute(
      range_km = range_m / 1000,
      method,
      metric = "RMSE",
      value = RMSE
    ),
  
  validation_scores %>%
    transmute(
      range_km = range_m / 1000,
      method,
      metric = "MAE",
      value = MAE
    )
)

range_plot <- ggplot(
  score_plot_data,
  aes(
    x = range_km,
    y = value,
    colour = method,
    group = method
  )
) +
  geom_vline(
    xintercept = range_m / 1000,
    linetype = "dotted",
    colour = "grey60"
  ) +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(
    breaks = range_values / 1000
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Sensitivity to spatial range",
    x = "Spatial range (km)",
    y = "Held-out error (kg C/m²)",
    colour = NULL,
    caption = "Dotted line: current range. Lower errors are better."
  ) +
  theme_minimal(base_size = 12)

print(range_plot)

ggsave(
  file.path(output_folder, "range_comparison.png"),
  range_plot,
  width = 9,
  height = 4.5,
  dpi = 300
)

cat(
  "\n5E complete. Current map range remains:",
  range_m, "metres.\n"
)



#==========================================
#----- 5F: Final summary and results -------
#==========================================

#----- Summarise the current run -----

score <- current_scores[
  current_scores$method == "Bayesian",
]

stopifnot(nrow(score) == 1)

# Stock means describe the valid mapped cells.

run_summary <- data.frame(
  target_depth_cm = target_depth_cm,
  range_m = range_m,
  observation_log_sd = observation_log_sd,
  
  n_cores = nrow(cores),
  n_plots = n_distinct(cores$plot_id),
  n_mapped_cells = nrow(posterior_values),
  
  prior_mean_stock = mean(posterior_values$prior_mean),
  updated_mean_stock = mean(posterior_values$posterior_mean),
  mean_update = mean(posterior_values$change_from_prior),
  
  heldout_RMSE = score$RMSE,
  heldout_MAE = score$MAE,
  heldout_bias = score$bias,
  coverage95_percent = score$coverage95_percent
)

write_csv(
  run_summary,
  file.path(output_folder, "run_summary.csv")
)


#----- Reproducible one-line summary -----

summary_line <- sprintf(
  paste0(
    "The 0-%d cm model used %d cores from %d plots ",
    "and a %.1f km spatial range. ",
    "Mean mapped stock was %.2f kg C/m2 ",
    "(prior: %.2f kg C/m2); ",
    "held-out RMSE was %.2f kg C/m2 ",
    "and 95%% predictive coverage was %.1f%%."
  ),
  target_depth_cm,
  nrow(cores),
  n_distinct(cores$plot_id),
  range_m / 1000,
  mean(posterior_values$posterior_mean),
  mean(posterior_values$prior_mean),
  score$RMSE,
  score$coverage95_percent
)

writeLines(
  c(
    summary_line,
    "",
    "Interpretation:",
    "Map means cover cells with usable prior information.",
    "Differences from the prior describe a model update.",
    "They do not measure carbon accumulation over time.",
    "Intervals are conditional on the recorded model assumptions.",
    "Prior SD was approximated from a 90% prediction-interval width.",
    "Pointwise intervals do not establish uncertainty in AOI totals."
  ),
  file.path(output_folder, "summary.txt")
)

cat("\n", summary_line, "\n")


#----- Record the input files -----

input_manifest <- data.frame(
  input = names(input_files),
  path = unname(input_files),
  md5 = unname(tools::md5sum(input_files))
)

write_csv(
  input_manifest,
  file.path(output_folder, "input_manifest.csv")
)


#----- Save a reusable results bundle -----

# Store ordinary R objects and grid geometry.
# The raster itself remains in bayesian_carbon_maps.tif.

results <- list(
  settings = model_settings,
  summary = run_summary,
  inputs = input_manifest,
  
  cores = cores,
  core_xy = core_xy,
  observation_var = observation_var,
  
  grid = list(
    nrow = nrow(bayesian_maps),
    ncol = ncol(bayesian_maps),
    extent = as.vector(ext(bayesian_maps)),
    crs = crs(bayesian_maps)
  ),
  
  grid_predictions = posterior_values,
  posterior_log_mean = posterior_log_mean,
  posterior_log_var = posterior_log_var,
  
  validation_predictions = validation_predictions,
  validation_scores = validation_scores
)

saveRDS(
  results,
  file.path(output_folder, "results.rds")
)


#----- Record software versions -----

capture.output(
  sessionInfo(),
  file = file.path(output_folder, "session_info.txt")
)

cat(
  "\nScript 5 complete. Results saved in:\n",
  normalizePath(output_folder),
  "\n"
)