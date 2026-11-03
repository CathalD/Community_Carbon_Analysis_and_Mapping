# Whole-core holdouts and sensitivity for the workshop raw-core update.
library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/workshop_gp.R")
folder <- "outputs/script5/sothe/stock_0_30/raw_cores"
state <- readRDS(file.path(folder, "raw_core_results.rds"))
cores <- state$cores
cv <- bind_rows(lapply(seq_len(nrow(cores)), function(i) {
  train <- setdiff(seq_len(nrow(cores)), i)
  p <- workshop_gp(state$core_xy[train,,drop=FALSE], cores$stock[train],
    cores$prior_mean[train], cores$prior_sd[train], state$core_xy[i,,drop=FALSE],
    cores$prior_mean[i], cores$prior_sd[i], fs_config$range_m, fs_config$observation_log_sd)
  data.frame(core_id = cores$core_id[i], observed = cores$stock[i], prior = p$prior_mean,
    predicted = p$posterior_mean, lower95 = p$lower95, upper95 = p$upper95,
    error = p$posterior_mean-cores$stock[i])
}))
write_csv(cv, file.path(folder, "leave_one_core_out.csv"))
write_csv(data.frame(n_cores = nrow(cv), prior_RMSE = sqrt(mean((cv$prior-cv$observed)^2)),
  updated_RMSE = sqrt(mean(cv$error^2)), MAE = mean(abs(cv$error)), bias = mean(cv$error),
  conditional_interval_coverage = mean(cv$observed >= cv$lower95 & cv$observed <= cv$upper95),
  status = "Eight individual-core holdouts; uncertainty is an S100-derived S30 sensitivity proxy"),
  file.path(folder, "validation_summary.csv"))
prior <- rast("outputs/script5/sothe/stock_0_30/prepared_prior.tif")
grid <- as.data.frame(prior, xy = TRUE, cells = TRUE, na.rm = TRUE)
scenarios <- expand.grid(range_m = c(1000, 2000, 4000), sd_multiplier = c(.5, 1, 2))
results <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  z <- scenarios[i,]
  p <- workshop_gp(state$core_xy, cores$stock, cores$prior_mean, cores$prior_sd*z$sd_multiplier,
    as.matrix(grid[c("x", "y")]), grid$prior_mean, grid$prior_sd*z$sd_multiplier,
    z$range_m, fs_config$observation_log_sd)
  if (z$range_m == fs_config$range_m && z$sd_multiplier == 1) {
    # Factored helper must reproduce the unchanged workshop implementation.
    original <- state$predictions
    stopifnot(identical(as.integer(original$cell), as.integer(grid$cell)),
      max(abs(p$posterior_mean-original$posterior_mean)) < 1e-8,
      max(abs(p$posterior_sd-original$posterior_sd)) < 1e-8)
  }
  data.frame(z, mean_updated = mean(p$posterior_mean), mean_sd = mean(p$posterior_sd))
}))
write_csv(results, file.path(folder, "range_and_uncertainty_sensitivity.csv"))
