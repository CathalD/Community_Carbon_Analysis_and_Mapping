# Workshop Script 3 adapted to the completed Fort Severn covariate exports.
# Original lm / leave-one-core-out / AICc / mapping steps are retained below.
library(terra)
library(sf)
library(readr)
library(dplyr)
library(ggplot2)
library(jsonlite)
source("config.R")
stock_variable <- "stock_0_30"
output_folder <- file.path("outputs/script3", stock_variable)
dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)
cores <- read_csv("outputs/script1/mapping_data.csv", show_col_types=FALSE) %>%
  rename(core_id=`Core ID`, longitude=Longitude, latitude=Latitude,
    stock_0_30=`C Stock 0-30cm`, stock_0_50=`C Stock 0-50cm`,
    stock_0_100=`C Stock 0-100cm`, stock_total=`Total Carbon Stock`)
cores$stock <- cores[[stock_variable]]
covariates <- rast("data/covariates/all_predictors_1km.tif")
aoi <- project(vect("data/boundary/aoi.geojson"), crs(covariates))
core_points <- project(vect(cores, geom=c("longitude","latitude"), crs="EPSG:4326"), crs(covariates))
# Fit to native core extraction, rather than replacing observations with cell-centre values.
core_values <- read_csv("data/covariates/native_core_covariates.csv", show_col_types=FALSE)
stopifnot(!anyDuplicated(core_values$core_id), setequal(cores$core_id, core_values$core_id))
model_data <- cores %>% select(core_id, stock) %>% left_join(core_values, by="core_id")
predictor_names <- setdiff(names(core_values), "core_id")
core_checks <- data.frame(predictor=predictor_names,
  missing_cores=vapply(model_data[predictor_names], function(x) sum(!is.finite(x)), integer(1)),
  distinct_values=vapply(model_data[predictor_names], function(x) length(unique(x[is.finite(x)])), integer(1)))
write_csv(model_data, file.path(output_folder,"core_covariates.csv"))
write_csv(core_checks, file.path(output_folder,"core_covariate_checks.csv"))
# An outcome-blind availability check retains all eight cores for fair model comparisons.
model_formulas <- list(); model_availability <- list()
for (model in names(fs_config$covariate_models)) {
  vars <- fs_config$covariate_models[[model]]
  supported <- all(vars %in% predictor_names) && all(vapply(vars, function(v)
    all(is.finite(model_data[[v]])) && length(unique(model_data[[v]])) >= 3, logical(1)))
  model_availability[[model]] <- data.frame(model=model, predictors=paste(vars,collapse=" + "),
    included=supported, rule="Every predictor finite at all eight cores; at least three distinct values")
  if (supported) model_formulas[[model]] <- reformulate(if(length(vars)) vars else "1", response="stock")
}
write_csv(bind_rows(model_availability), file.path(output_folder,"candidate_availability.csv"))
stopifnot(length(model_formulas) > 0)
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


