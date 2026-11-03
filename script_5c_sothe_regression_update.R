# Fourth Sothe branch: selected Script 3 regression map + workshop prior fusion.
# The pixelwise log-Gaussian likelihood and moment calculations are the same
# equations as Script 5's area-based and IDW branches. The working error is
# the selected model's whole-core holdout RMSE, not nested selection error.
library(terra)
library(readr)
library(ggplot2)

out <- "outputs/script5/sothe/stock_0_30/regression"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
prior_path <- "outputs/script5/sothe/stock_0_30/prepared_prior.tif"
regression_path <- "outputs/script3/stock_0_30/selected_carbon_stock.tif"
selection_path <- "outputs/script3/stock_0_30/selected_model.csv"
cv_path <- "outputs/script3/stock_0_30/leave_one_core_out.csv"
fields_path <- "outputs/script1/mapping_data.csv"
stopifnot(all(file.exists(c(prior_path, regression_path, selection_path, cv_path, fields_path))))

prior <- rast(prior_path)
stopifnot(all(c("prior_mean", "prior_sd") %in% names(prior)))
selected <- read_csv(selection_path, show_col_types = FALSE)
stopifnot(nrow(selected) == 1, selected$stock_variable == "stock_0_30")
model_name <- selected$model[[1]]
holdouts <- read_csv(cv_path, show_col_types = FALSE)
holdouts <- holdouts[holdouts$model == model_name, ]
fields <- read_csv(fields_path, show_col_types = FALSE)
stopifnot(nrow(holdouts) == nrow(fields), !anyDuplicated(holdouts$core_id),
  setequal(holdouts$core_id, fields[["Core ID"]]),
  all(is.finite(holdouts$observed)), all(is.finite(holdouts$predicted)),
  all(holdouts$observed > 0))
holdouts$error <- holdouts$predicted - holdouts$observed
error_sd <- sqrt(mean(holdouts$error^2))
stopifnot(is.finite(error_sd), error_sd > 0,
  abs(error_sd - selected$RMSE[[1]]) < 1e-8)

regression <- rast(regression_path)
if (!compareGeom(regression, prior[[1]], stopOnError = FALSE)) {
  regression <- project(regression, prior[[1]], method = "bilinear")
}
regression <- mask(regression, prior$prior_mean)
names(regression) <- "regression_stock"
stopifnot(compareGeom(regression, prior[[1]], stopOnError = FALSE))
writeRaster(regression, file.path(out, "regression_stock.tif"), overwrite = TRUE)

result <- as.data.frame(c(prior, regression), cells = TRUE, xy = TRUE, na.rm = TRUE)
eligible <- is.finite(result$prior_mean) & result$prior_mean > 0 &
  is.finite(result$prior_sd) & result$prior_sd > 0 &
  is.finite(result$regression_stock) & result$regression_stock > 0
result <- result[eligible, ]
stopifnot(nrow(result) > 0)

# Workshop Script 5 pointwise log-Gaussian update. Its working likelihood
# treats log(mapped stock) as an unbiased observation of the latent log stock.
prior_var <- log1p((result$prior_sd / result$prior_mean)^2)
prior_mu <- log(result$prior_mean) - prior_var / 2
likelihood_var <- log1p((error_sd / result$regression_stock)^2)
weight <- prior_var / (prior_var + likelihood_var)
post_mu <- prior_mu + weight * (log(result$regression_stock) - prior_mu)
post_var <- (1 - weight) * prior_var
result$posterior_log_mean <- post_mu
result$posterior_log_var <- post_var
result$data_weight <- weight
result$posterior_mean <- exp(post_mu + post_var/2)
result$posterior_sd <- sqrt(expm1(post_var) * exp(2*post_mu + post_var))
result$lower95 <- exp(post_mu - qnorm(.975)*sqrt(post_var))
result$upper95 <- exp(post_mu + qnorm(.975)*sqrt(post_var))
result$change_from_prior <- result$posterior_mean - result$prior_mean
result$log_variance_reduction <- 1-post_var/prior_var

# Independently verify the precision-weighted form and lognormal moments.
precision_var <- 1 / (1/prior_var + 1/likelihood_var)
precision_mu <- precision_var * (prior_mu/prior_var + log(result$regression_stock)/likelihood_var)
checks <- data.frame(
  check = c("Every field core has a selected-model holdout",
    "Posterior log variance nonnegative and no larger than prior",
    "Precision and weight forms give identical posterior",
    "Moment reconstruction matches mapped posterior mean",
    "Finite positive posterior means and SDs",
    "Intervals bracket posterior median"),
  passed = c(nrow(holdouts) == 8,
    all(post_var >= 0 & post_var <= prior_var),
    max(abs(post_var-precision_var)) < 1e-12 && max(abs(post_mu-precision_mu)) < 1e-12,
    max(abs(exp(precision_mu+precision_var/2)-result$posterior_mean)) < 1e-10,
    all(is.finite(as.matrix(result))) && all(result$posterior_mean > 0 & result$posterior_sd > 0),
    all(result$lower95 <= exp(post_mu) & exp(post_mu) <= result$upper95)))
write_csv(checks, file.path(out, "verification.csv"))
stopifnot(all(checks$passed))

map_columns <- c("prior_mean", "prior_sd", "regression_stock", "posterior_mean", "posterior_sd",
  "lower95", "upper95", "change_from_prior", "log_variance_reduction", "data_weight")
maps <- rast(prior[[1]], nlyrs = length(map_columns))
names(maps) <- map_columns
v <- matrix(NA_real_, ncell(maps), length(map_columns))
v[result$cell, ] <- as.matrix(result[map_columns])
values(maps) <- v
writeRaster(maps, file.path(out, "posterior_maps.tif"), overwrite = TRUE)
write_csv(result, file.path(out, "predictions.csv"))
write_csv(holdouts, file.path(out, "error_checks.csv"))
write_csv(selected, file.path(out, "selected_regression_model.csv"))
write_csv(data.frame(method = "Covariate-regression update", selected_model = model_name,
  n_cores = nrow(holdouts), prediction_cells = nrow(result),
  mean_prior_stock_kg_m2 = mean(result$prior_mean),
  mean_updated_stock_kg_m2 = mean(result$posterior_mean),
  working_error_sd_kg_m2 = error_sd,
  error_source = "Whole-core holdout RMSE; selected formula fixed, not nested model-selection validation",
  model = "Approximate pointwise log-Gaussian fusion; no joint covariance",
  likelihood = "Log-unbiased error; plug-in log variance",
  regression_map = regression_path,
  uncertainty_status = "S100-derived absolute S30 SD proxy; conditional sensitivity intervals"),
  file.path(out, "summary.csv"))

p <- ggplot(result, aes(x, y, fill = posterior_mean)) + geom_raster() +
  scale_fill_viridis_c(name = "kg C/m²") + coord_equal() + theme_minimal() +
  labs(title = "Sothe 0–30 cm: regression-map prior update",
    subtitle = paste("Selected Script 3 model:", model_name), x = "Easting (m)", y = "Northing (m)")
ggsave(file.path(out, "posterior_mean.png"), p, width = 8, height = 7, dpi = 180)
writeLines(c(
  "The regression branch uses the selected Script 3 map and the same pointwise fusion equations as workshop Script 5.",
  paste("Selected model:", model_name),
  "A Mean_only selection remains an intercept-only regression; an environmental effect is not forced into the model.",
  "The holdout RMSE fixes the selected formula; it does not include nested model selection.",
  "This and the other three updates reuse the same eight cores and are alternative estimates, not independent evidence.",
  "The uncertainty is conditional on the cross-depth S100-to-S30 uncertainty proxy and working likelihood."),
  file.path(out, "interpretation.txt"))
message("Sothe regression prior update completed: ", model_name)
