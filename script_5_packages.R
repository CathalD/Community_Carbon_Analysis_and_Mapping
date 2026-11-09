# Optional numerical comparison after Script 5. Not part of the workshop run.
# Original full package experiment: archive/development/script_5_packages.R.
library(terra)
library(readr)
library(gstat)
library(mvtnorm)
library(sp)
source("config.R")
terra::projNetwork(FALSE)

folder <- file.path("outputs/script5", paste0("stock_0_", bayesian_depth_cm))
baseline <- readRDS(file.path(folder, "raw_cores/model.rds"))
written <- read_csv(file.path(folder, "raw_cores/predictions.csv"), show_col_types = FALSE)
prior <- rast(file.path("outputs/script4", paste0("prior_moments_0_", bayesian_depth_cm, "cm.tif")))
grid <- values(prior, mat = TRUE)[written$cell, , drop = FALSE]
core_var <- log1p((baseline$cores$prior_sd / baseline$cores$prior_mean)^2)
core_mean <- log(baseline$cores$prior_mean) - core_var / 2
projection <- sp::CRS(SRS_string = crs(prior))
training <- sp::SpatialPointsDataFrame(
  baseline$core_xy,
  data = data.frame(
    residual = (log(baseline$cores$stock) - core_mean) / sqrt(core_var),
    noise_var = baseline$observation_log_sd^2 / core_var
  ), proj4string = projection
)
targets <- sp::SpatialPointsDataFrame(
  as.matrix(written[, c("x", "y")]),
  data = data.frame(noise_var = rep(0, nrow(written))), proj4string = projection
)
raw_covariance <- function(x, y = x) {
  C <- exp(-sp::spDists(sp::coordinates(x), sp::coordinates(y)) / baseline$range_m)
  if (identical(x, y)) diag(C) <- diag(C) + x$noise_var
  C
}
package <- gstat::krige0(
  residual ~ 1, data = training, newdata = targets, model = raw_covariance,
  beta = 0, computeVar = TRUE, fullCovariance = FALSE
)
comparisons <- list(data.frame(
  method = "Raw cores: gstat",
  max_log_mean_difference = max(abs(
    grid[, "log_mean"] + sqrt(grid[, "log_variance"]) * as.vector(package$pred) -
      written$posterior_log_mean
  )),
  max_log_variance_difference = max(abs(
    grid[, "log_variance"] * as.vector(package$var) - written$posterior_log_var
  ))
))

for (method in c("area_based", "idw")) {
  written <- read_csv(file.path(folder, method, "predictions.csv"), show_col_types = FALSE)
  grid <- values(prior, mat = TRUE)[written$cell, , drop = FALSE]
  local <- rast(file.path("outputs/script2", paste0(method, "_log_moments_0_", bayesian_depth_cm, "cm.tif")))
  local_values <- values(local, mat = TRUE)[written$cell, , drop = FALSE]
  q <- local_values[, "log_variance"] / grid[, "log_variance"]
  covariance <- mvtnorm::syMatrices(rbind(1 + q, rep(1, length(q)), rep(1, length(q))), diag = TRUE)
  conditional <- mvtnorm::cond_mvnorm(
    chol = chol(covariance), which_given = 1L,
    given = matrix((local_values[, "log_mean"] - grid[, "log_mean"]) /
                     sqrt(grid[, "log_variance"]), nrow = 1L)
  )
  comparisons[[method]] <- data.frame(
    method = paste(method, "mvtnorm"),
    max_log_mean_difference = max(abs(
      grid[, "log_mean"] + sqrt(grid[, "log_variance"]) * as.numeric(conditional$mean) -
        written$posterior_log_mean
    )),
    max_log_variance_difference = max(abs(
      grid[, "log_variance"] * as.numeric(mvtnorm::Tcrossprod(conditional$chol, diag_only = TRUE)) -
        written$posterior_log_var
    ))
  )
}
cat("\nNumerical implementation differences; these do not assess predictive accuracy:\n")
print(do.call(rbind, comparisons))
