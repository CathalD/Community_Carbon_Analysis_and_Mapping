library(terra)
library(sf)
library(readr)
library(dplyr)
library(spBayes)

# 1. Load Data & Prior Rasters
target_depth_cm <- bayesian_depth_cm
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")

shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask

# Load National Prior Raster
prior <- rast(file.path("outputs/script4", paste0("prior_moments_0_", target_depth_cm, "cm.tif")))
prior <- mask(prior, analysis_mask)
prior_values <- values(prior, mat = TRUE)
eligible_cells <- which(apply(is.finite(prior_values), 1, all) & prior_values[, "prior_mean"] > 0)

grid_xy <- xyFromCell(template, eligible_cells)
grid_mean <- prior_values[eligible_cells, "prior_mean"]

# Load Local Field Cores
raw_cores <- read_csv(mapping_file, show_col_types = FALSE)
cores <- data.frame(
  core_id = raw_cores[["Core ID"]], 
  longitude = raw_cores$Longitude, 
  latitude = raw_cores$Latitude,
  stock = raw_cores[[stock_column]]
)
cores <- cores[is.finite(cores$stock) & cores$stock > 0, ]
core_points <- project(vect(cores, geom = c("longitude", "latitude"), crs = coordinate_crs), crs(template))
core_xy <- crds(core_points)

# Extract National Prior at Core Locations
depth_key <- as.character(target_depth_cm)
native_prior <- rast(prior_mean_files[[depth_key]])
cores$log_prior_mean <- log(extract(native_prior, project(core_points, crs(native_prior)), ID = FALSE)[, 1])

# Format Coordinates as Matrices for spBayes
coords      <- as.matrix(core_xy)
grid_coords <- as.matrix(grid_xy)

# 2. Define Priors (Anchored to National Expectations)
priors <- list(
  beta.normal   = list(mean = c(0, 1), var = c(10, 10)), # Intercept 0, Slope 1 with national prior
  sigma.sq.IG   = c(2, 1),                               # Spatial process variance (sill)
  tau.sq.IG     = c(2, 1),                               # Nugget error variance
  phi.Unif      = c(0.0001, 0.01)                        # Spatial decay range bounds
)

# MCMC Starting Values & Tuning parameters
starting <- list(phi = 0.001, sigma.sq = 1, tau.sq = 0.5, beta = c(0, 1))
tuning   <- list(phi = 0.001, sigma.sq = 0.05, tau.sq = 0.05, beta = c(0.01, 0.01))

# 3. Fit the Hierarchical Spatial Linear Model (spLM)
set.seed(42)
cat("Fitting Gaussian Process via spBayes MCMC...\n")
sp_model <- spLM(
  formula = log(stock) ~ log_prior_mean,
  data = cores,
  coords = coords,
  cov.model = "exponential",
  priors = priors,
  starting = starting,
  tuning = tuning,
  n.samples = 10000,
  verbose = FALSE
)

# 4. Recover Spatial Random Effects (Drop 5,000 burn-in iterations)
sp_model_recovered <- spRecover(sp_model, start = 5001, thin = 2, verbose = FALSE)

# 5. Predict Across the Shared Raster Grid
cat("Generating posterior predictive distributions across raster grid...\n")
pred_covars <- cbind(intercept = 1, log_prior_mean = log(grid_mean))
pred_results <- spPredict(
  sp_model_recovered,
  pred.coords = grid_coords,
  pred.covars = pred_covars,
  start = 5001,
  thin = 2
)

# Extract MCMC Predictive Samples (Rows = Grid Pixels, Cols = MCMC Iterations)
pred_samples <- pred_results$p.y.predictive.samples
sp_log_mean <- apply(pred_samples, 1, mean)
sp_log_var  <- apply(pred_samples, 1, var)

# 6. Export Back to Rasters (With Unbiased Log-Normal Correction)
results <- data.frame(
  cell = eligible_cells, x = grid_xy[, 1], y = grid_xy[, 2],
  posterior_mean = exp(sp_log_mean + sp_log_var / 2),
  posterior_sd   = sqrt(expm1(sp_log_var) * exp(2 * sp_log_mean + sp_log_var))
)

maps <- rast(template, nlyrs = 2)
names(maps) <- c("gp_posterior_mean", "gp_posterior_sd")
cells <- matrix(NA_real_, ncell(template), 2)
cells[eligible_cells, ] <- as.matrix(results[, c("posterior_mean", "posterior_sd")])
values(maps) <- cells
maps <- mask(maps, analysis_mask)

writeRaster(maps[[1]], "outputs/gp_posterior_mean.tif", overwrite = TRUE)
writeRaster(maps[[2]], "outputs/gp_posterior_sd.tif", overwrite = TRUE)
cat("\nGaussian Process mapping complete and saved to raster!\n")


plot(maps)
