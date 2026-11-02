# ==============================================

# ------- Part  3- Spatial Mapping Continued ---------

# ===========================================

# Install missing packages once:
# install.packages(c("terra", "sf", "readr", "dplyr", "ggplot2", "jsonlite"))

library(terra)
library(sf)
library(readr)
library(dplyr)
library(ggplot2)
library(jsonlite)

#---------
# Step 1 - Load data and covaraites
# -------

# Input files
mapping_file <- "outputs/script1/mapping_data.csv"
aoi_file <- "data/boundary/aoi.geojson"

elevation_file <- "data/covariates/topography.tif"
water_file <- "data/covariates/water.tif"

# To start we will run a 2-D model
# Starting with 0-30cm stock

stock_variable <- "stock_0_30"

# Create output folder
output_folder <- file.path("outputs/script3", stock_variable)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Read the field stocks from Script 1
cores <- read_csv(mapping_file, show_col_types = FALSE) %>%
  rename(
    core_id = `Core ID`,
    longitude = Longitude,
    latitude = Latitude,
    stock_0_30 = `C Stock 0-30cm`,
    stock_0_50 = `C Stock 0-50cm`,
    stock_0_100 = `C Stock 0-100cm`,
    stock_total = `Total Carbon Stock`
  )

cores$stock <- cores[[stock_variable]]

# Read the two rasters
elevation <- rast(elevation_file)
water <- rast(water_file)


# GIS housekeeping - check names and CRS

names(elevation) <- "elevation_m"
names(water) <- "water_occurrence_pct"

if (!compareGeom(elevation, water, stopOnError = FALSE)) {
  stop("Elevation and water do not share the same raster grid.")
}

# Display dimensions, resolution and coordinate system
print(elevation)
print(water)

# Display the range of values
global(c(elevation, water), c("min", "max"), na.rm = TRUE)

# Read the boundary and match it to the raster coordinate system
aoi <- vect(aoi_file)
aoi <- project(aoi, crs(elevation))

# Convert field coordinates into map points
core_points <- vect(
  cores[, c("core_id", "longitude", "latitude", "stock")],
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

core_points <- project(core_points, crs(elevation))

# Check the selected field stocks
print(cores %>% select(core_id, stock))

# Display both maps with the boundary and core locations
old_par <- par(mfrow = c(1, 2))

plot(elevation, main = "Elevation (m)")
plot(aoi, add = TRUE, border = "black")
points(core_points, pch = 21, bg = "red", cex = 1)

plot(water, main = "Water occurrence (%)")
plot(aoi, add = TRUE, border = "black")
points(core_points, pch = 21, bg = "red", cex = 1)

par(old_par)


#===============================================================
# Step 3C - Prep the Covariate Files
#===============================================================

# Terrain derivatives from the elevation raster
terrain_layers <- terrain(
  elevation,
  v = c("slope", "aspect", "TPI", "TRI"),
  unit = "degrees",
  neighbors = 8
)

names(terrain_layers) <- c(
  "slope_deg",
  "aspect_deg",
  "tpi_m",
  "tri_m"
)

# Express aspect as northness and eastness for modelling.
# Flat cells receive zero for both direction variables.
aspect_radians <- terrain_layers$aspect_deg * pi / 180

northness <- ifel(
  terrain_layers$slope_deg == 0,
  0,
  cos(aspect_radians)
)

eastness <- ifel(
  terrain_layers$slope_deg == 0,
  0,
  sin(aspect_radians)
)

names(northness) <- "northness"
names(eastness) <- "eastness"

# Distance to water
water_threshold <- 50

water_cells <- ifel(water >= water_threshold, 1, NA)

if (global(water_cells, "notNA")[1, 1] == 0) {
  stop("No water cells meet the selected occurrence threshold.")
}

distance_to_water <- distance(water_cells, unit = "m")

# Retain the original water raster's data coverage
distance_to_water <- mask(distance_to_water, water)
names(distance_to_water) <- "distance_to_water_m"

# Combine the two original layers and their derivatives
covariates <- c(
  elevation,
  terrain_layers,
  northness,
  eastness,
  water,
  distance_to_water
)

# Apply the study boundary after deriving the predictors
covariates <- mask(covariates, aoi, touches = FALSE)

# Inspect the resulting layers
print(covariates)
global(covariates, c("min", "max"), na.rm = TRUE)

plot(covariates)

# Save the shared predictor layers
covariate_output <- "outputs/script3/covariates"
dir.create(covariate_output, recursive = TRUE, showWarnings = FALSE)

writeRaster(
  covariates,
  file.path(covariate_output, "terrain_water.tif"),
  overwrite = TRUE
)


# ---------------
# Extract Covariate Values at Field Locations
# --------------------

core_values <- terra::extract(
  covariates,
  core_points,
  ID = FALSE
)

# One row per core, containing stock and environmental predictors
model_data <- bind_cols(
  cores %>% select(core_id, stock),
  core_values
)

print(model_data, width = Inf)

# Check missing values and variation between cores
core_checks <- data.frame(
  predictor = names(core_values),
  missing_cores = colSums(is.na(core_values)),
  distinct_values = sapply(
    core_values,
    function(x) length(unique(x[!is.na(x)]))
  )
)

print(core_checks)

# Save the extracted data and checks
write_csv(
  model_data,
  file.path(output_folder, "core_covariates.csv")
)

write_csv(
  core_checks,
  file.path(output_folder, "core_covariate_checks.csv")
)

# Compare each predictor visually with field carbon stock
library(ggplot2)

relationship_data <- bind_rows(lapply(
  names(core_values),
  function(predictor) {
    data.frame(
      core_id = model_data$core_id,
      stock = model_data$stock,
      predictor = predictor,
      value = model_data[[predictor]]
    )
  }
))

relationship_plot <- ggplot(
  relationship_data,
  aes(x = value, y = stock)
) +
  geom_point(size = 2.5, na.rm = TRUE) +
  facet_wrap(~ predictor, scales = "free_x", ncol = 3) +
  labs(
    title = "Field carbon stocks and environmental predictors",
    x = "Predictor value",
    y = "Carbon stock (kg C/m²)"
  ) +
  theme_minimal()

print(relationship_plot)

# Remember this data is randomly generated - So it shoudl look random

ggsave(
  file.path(output_folder, "predictor_relationships.png"),
  relationship_plot,
  width = 11, height = 8, dpi = 200
)


#===============================================================
# Step 3E — Use this information to build a set of simple models we can use to
# help extrapolate the values across the AOI with some spatial knowledge
#===============================================================

model_formulas <- list(
  Mean_only = stock ~ 1,
  Elevation = stock ~ elevation_m,
  Slope = stock ~ slope_deg,
  TPI = stock ~ tpi_m,
  Water_distance = stock ~ distance_to_water_m,
  Combined = stock ~ elevation_m + distance_to_water_m
)

# Use the same complete observations for every model
required_predictors <- unique(unlist(lapply(
  model_formulas,
  function(f) all.vars(f)[-1]
)))

required_columns <- c("stock", required_predictors)

keep <- Reduce(
  `&`,
  lapply(model_data[required_columns], is.finite)
)

training_data <- model_data[keep, ]
excluded_data <- model_data[!keep, ]

# Allow enough observations for fitting and held-out predictions
max_predictors <- max(vapply(
  model_formulas,
  function(f) length(all.vars(f)) - 1L,
  integer(1)
))

if (nrow(training_data) < max_predictors + 3) {
  stop("Too few complete cores for these models and cross-validation.")
}

cat("Cores used by every model:", nrow(training_data), "\n")
cat("Cores excluded:", nrow(excluded_data), "\n")

write_csv(training_data, file.path(output_folder, "training_data.csv"))
write_csv(excluded_data, file.path(output_folder, "excluded_cores.csv"))


#===============================================================
# Step 3F — Fit the models and Run Unceratinty analysis to evaluate effective of the model
#===============================================================

leave_one_core_out <- function(formula, data) {
  predicted <- numeric(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    fit <- lm(
      formula,
      data = data[-i, ],
      na.action = na.fail
    )
    
    if (anyNA(coef(fit))) {
      stop(
        "A held-out fit has constant or dependent predictors: ",
        paste(deparse(formula), collapse = " ")
      )
    }
    
    predicted[i] <- as.numeric(predict(
      fit,
      newdata = data[i, , drop = FALSE]
    ))
  }
  
  data.frame(
    core_id = data$core_id,
    observed = data$stock,
    predicted = predicted
  )
}

fitted_models <- list()
validation_results <- list()
model_scores <- list()

for (model_name in names(model_formulas)) {
  formula <- model_formulas[[model_name]]
  
  # Fit the final model using all eligible cores
  fitted_models[[model_name]] <- lm(
    formula,
    data = training_data,
    na.action = na.fail
  )
  
  # Generate held-out predictions
  validation <- leave_one_core_out(formula, training_data)
  validation$model <- model_name
  
  # Positive error means overprediction
  validation$error <- validation$predicted - validation$observed
  
  total_variation <- sum(
    (validation$observed - mean(validation$observed))^2
  )
  
  model_scores[[model_name]] <- data.frame(
    model = model_name,
    n_cores = nrow(validation),
    RMSE = sqrt(mean(validation$error^2)),
    MAE = mean(abs(validation$error)),
    bias = mean(validation$error),
    R2 = if (total_variation > 0) {
      1 - sum(validation$error^2) / total_variation
    } else {
      NA_real_
    }
  )
  
  validation_results[[model_name]] <- validation
}

# Combine and display the results
model_statistics <- bind_rows(model_scores) %>%
  arrange(RMSE)

validation_predictions <- bind_rows(validation_results)

print(model_statistics)

write_csv(
  model_statistics,
  file.path(output_folder, "model_statistics.csv")
)

write_csv(
  validation_predictions,
  file.path(output_folder, "leave_one_core_out.csv")
)

saveRDS(
  fitted_models,
  file.path(output_folder, "fitted_models.rds")
)

# Observed versus held-out predicted stocks
validation_plot <- ggplot(
  validation_predictions,
  aes(x = observed, y = predicted)
) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", colour = "grey50"
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ model, ncol = 3) +
  coord_equal() +
  labs(
    title = "Observed and held-out predicted carbon stocks",
    x = "Observed stock (kg C/m²)",
    y = "Predicted stock (kg C/m²)"
  ) +
  theme_minimal()

print(validation_plot)

ggsave(
  file.path(output_folder, "observed_vs_predicted.png"),
  validation_plot,
  width = 10, height = 7, dpi = 200
)

#===============================================================
# Step 3G - Apply the models back into spatial outputs
#===============================================================

prediction_maps <- list()
map_summaries <- list()

for (model_name in names(fitted_models)) {
  cat("Mapping:", model_name, "\n")
  
  fit <- fitted_models[[model_name]]
  predictors <- all.vars(formula(fit))[-1]
  
  model_folder <- file.path(output_folder, model_name)
  dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)
  
  if (length(predictors) == 0) {
    
    # Mean-only baseline: constant stock within the mapped area
    prediction <- ifel(
      !is.na(covariates$elevation_m),
      as.numeric(coef(fit)[1]),
      NA
    )
    
    outside_range <- ifel(!is.na(prediction), 0, NA)
    
  } else {
    
    # Supply only the layers used by this model
    model_rasters <- covariates[[predictors]]
    
    prediction <- terra::predict(
      model_rasters,
      fit,
      fun = stats::predict,
      na.rm = TRUE
    )
    
    # Identify predictor values beyond those observed at the cores
    first_predictor <- predictors[1]
    
    outside_range <-
      model_rasters[[first_predictor]] <
      min(training_data[[first_predictor]]) |
      model_rasters[[first_predictor]] >
      max(training_data[[first_predictor]])
    
    for (predictor in predictors[-1]) {
      outside_range <- outside_range |
        model_rasters[[predictor]] <
        min(training_data[[predictor]]) |
        model_rasters[[predictor]] >
        max(training_data[[predictor]])
    }
  }
  
  names(prediction) <- "carbon_stock_kg_m2"
  
  writeRaster(
    prediction,
    file.path(model_folder, "carbon_stock_kg_m2.tif"),
    overwrite = TRUE
  )
  
  # Keep the original predictions and save checks separately
  prediction_flags <- c(outside_range, prediction < 0)
  
  names(prediction_flags) <- c(
    "outside_training_range",
    "negative_prediction"
  )
  
  writeRaster(
    prediction_flags,
    file.path(model_folder, "prediction_flags.tif"),
    overwrite = TRUE
  )
  
  # Summarize each map
  map_stats <- global(
    prediction,
    c("min", "max", "notNA"),
    na.rm = TRUE
  )
  
  flag_counts <- global(
    prediction_flags,
    "sum",
    na.rm = TRUE
  )
  
  map_summaries[[model_name]] <- data.frame(
    model = model_name,
    minimum_stock = map_stats$min,
    maximum_stock = map_stats$max,
    valid_cells = map_stats$notNA,
    outside_training_range_cells = flag_counts$sum[1],
    negative_prediction_cells = flag_counts$sum[2]
  )
  
  prediction_maps[[model_name]] <- prediction
}

map_summary <- bind_rows(map_summaries)

print(map_summary)

write_csv(
  map_summary,
  file.path(output_folder, "map_summary.csv")
)


#---------------------------------------------------------------
# Display and save maps using one shared colour scale
#---------------------------------------------------------------

map_limits <- range(
  map_summary$minimum_stock,
  map_summary$maximum_stock,
  na.rm = TRUE
)

draw_model_map <- function(prediction, model_name) {
  plot(
    prediction,
    main = paste(model_name, stock_variable, sep = "\n"),
    range = map_limits,
    plg = list(title = "kg C/m²")
  )
  
  plot(aoi, add = TRUE, border = "black")
  points(core_points, pch = 21, bg = "white", cex = 0.9)
}

for (model_name in names(prediction_maps)) {
  prediction <- prediction_maps[[model_name]]
  model_folder <- file.path(output_folder, model_name)
  
  png(
    file.path(model_folder, "carbon_stock_map.png"),
    width = 1400, height = 1100, res = 150
  )
  
  draw_model_map(prediction, model_name)
  dev.off()
  
  # Also display in the R plotting window
  draw_model_map(prediction, model_name)
}


#===============================================================
# Step 3H — Evaluate and select the model for Script 4
#===============================================================

aic_results <- bind_rows(lapply(
  names(fitted_models),
  function(model_name) {
    fit <- fitted_models[[model_name]]
    
    n <- nobs(fit)
    K <- attr(logLik(fit), "df")  # Includes residual variance
    aic <- AIC(fit)
    
    aicc <- if (n > K + 1 && is.finite(aic)) {
      aic + (2 * K * (K + 1)) / (n - K - 1)
    } else {
      NA_real_
    }
    
    data.frame(
      model = model_name,
      parameters = K,
      AIC = aic,
      AICc = aicc
    )
  }
))

if (!any(is.finite(aic_results$AICc))) {
  stop("AICc could not be calculated for these models.")
}

# Combine information criteria with held-out prediction scores
model_evaluation <- model_statistics %>%
  left_join(aic_results, by = "model") %>%
  mutate(
    delta_AICc = AICc - min(AICc, na.rm = TRUE)
  ) %>%
  arrange(AICc, RMSE)

print(model_evaluation)

write_csv(
  model_evaluation,
  file.path(output_folder, "model_evaluation.csv")
)

# Display the relative AICc support
evaluation_plot <- model_evaluation %>%
  filter(is.finite(delta_AICc)) %>%
  ggplot(aes(
    x = delta_AICc,
    y = reorder(model, -delta_AICc)
  )) +
  geom_point(size = 3, colour = "#34778c") +
  labs(
    title = "Model evaluation",
    subtitle = "Lower ΔAICc indicates greater support",
    x = "ΔAICc",
    y = NULL
  ) +
  theme_minimal()

print(evaluation_plot)

ggsave(
  file.path(output_folder, "model_evaluation.png"),
  evaluation_plot,
  width = 7, height = 4, dpi = 200
)

# Select the lowest-AICc model
selected_record <- model_evaluation %>%
  filter(is.finite(AICc)) %>%
  slice(1) %>%
  mutate(
    stock_variable = stock_variable,
    selection_rule = "Lowest AICc"
  )

selected_model <- selected_record$model[1]

cat("Selected model:", selected_model, "\n")

# Locate the map already created in 3G
selected_map_file <- file.path(
  output_folder,
  selected_model,
  "carbon_stock_kg_m2.tif"
)

if (!file.exists(selected_map_file)) {
  stop("Run 3G before exporting the selected model map.")
}

# Save the fitted model, evaluation record and map
saveRDS(
  fitted_models[[selected_model]],
  file.path(output_folder, "selected_model.rds")
)

write_csv(
  selected_record,
  file.path(output_folder, "selected_model.csv")
)

file.copy(
  selected_map_file,
  file.path(output_folder, "selected_carbon_stock.tif"),
  overwrite = TRUE
)


