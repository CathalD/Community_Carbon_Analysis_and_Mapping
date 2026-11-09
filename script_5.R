# ===================================
# Part 5: Updating existing maps with local field data
# ===================================

library(terra)
library(readr)
library(mapview)

source("helpers/grid_helper.R")

# ----------------
# Step 1 - Set up the analysis
# ----------------

# Start with 100 cm. Run Scripts 2 and 4 first for the chosen depth.
target_depth_cm <- 100
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
output_folder <- file.path("outputs/script5", paste0("stock_0_", target_depth_cm))
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
source("config.R")

# Prediction points are the shared grid-cell centres inside the study area.
source("helpers/grid_helper.R") # Where we define our gis setting and predictions
shared <- load_analysis_grid()
analysis_mask <- shared$mask
predicted_points <- as.points(analysis_mask)
mapview(predicted_points)

# Import Field stocks calculated in Script 
cores <- read_csv(mapping_file)
field_stock <- cores[[stock_column]]
core_points <- vect(cores, geom = c("Longitude", "Latitude"), crs = coordinate_crs)
core_points <- project(core_points, crs(analysis_mask))

# Impor Prior from script 4: McMaster mean, SD and log transformed data
prior <- rast(file.path("outputs/script4", paste0("prior_moments_0_", target_depth_cm, "cm.tif")))
prior <- mask(prior, analysis_mask)
prior_mean <- prior$prior_mean
prior_sd <- prior$prior_sd

# Local predictions from Script 2: equal weights and distance weights.
area_based <- rast(file.path("outputs/script2", paste0("area_based_log_moments_0_", target_depth_cm, "cm.tif")))
idw <- rast(file.path("outputs/script2", paste0("idw_log_moments_0_", target_depth_cm, "cm.tif")))

# Lets view our input data

# Work on log carbon, where these models use Normal distributions.
# All v variables below are variances, not standard deviations.
m0 <- prior$log_mean
v0 <- prior$log_variance
m_area <- area_based$log_mean
v_area <- area_based$log_variance
m_idw <- idw$log_mean
v_idw <- idw$log_variance

mapview(m0) + mapview(area_based) + mapview(idw)


# ----------------
# Step 2 - Conjugate updating
# ----------------

# Bayes' theorem:
# p(Y | local data) = p(local data | Y) * p(Y) / p(local data)
#
# Y = log carbon stock.
# Prior: Y ~ Normal(m0, v0).
# Local estimate: yhat | Y ~ Normal(Y, v_local).
#
# Posterior variance = v0 * v_local / (v0 + v_local)
# Posterior mean = m0 + v0/(v0 + v_local) * (yhat - m0)
# The estimate with smaller variance receives more weight.
#
# Here the local map's prediction-error variance is a working likelihood variance:
# this is approximate map fusion, assuming unbiased, independent information.
# Area-based and IDW use the same cores; update the original prior separately.


#post_log_v =
#post_log_m = (m0 * v0``)


# A. Update with the area-based map.
post_log_variance <- v0 * v_area / (v0 + v_area)
post_log_mean <- m0 + v0 / (v0 + v_area) * (m_area - m0)

# Convert the posterior back to carbon units (kg C/m2).
# exp(log mean) is the median; the arithmetic mean includes variance/2.
post_stock <- exp(post_log_mean + post_log_variance / 2)
post_stock_sd <- sqrt((exp(post_log_variance) - 1) *
                        exp(2 * post_log_mean + post_log_variance))

mapview(post_stock)

# B. Update the same original prior with the IDW map.
post_log_variance_idw <- v0 * v_idw / (v0 + v_idw)
post_log_mean_idw <- m0 + v0 / (v0 + v_idw) * (m_idw - m0)
post_stock_idw <- exp(post_log_mean_idw + post_log_variance_idw / 2)
post_stock_sd_idw <- sqrt((exp(post_log_variance_idw) - 1) *
                            exp(2 * post_log_mean_idw + post_log_variance_idw))

# ----------------
# Step 3 - Show a distribution example
# ----------------

# Use the first cell with a prior value, and the area-based update.
example_cell <- which(!is.na(values(m0)))[1]
mu <- c(values(m0)[example_cell], values(m_area)[example_cell],
        values(post_log_mean)[example_cell])
sigma <- sqrt(c(values(v0)[example_cell], values(v_area)[example_cell],
                values(post_log_variance)[example_cell]))
x <- seq(min(mu - 4 * sigma), max(mu + 4 * sigma), length.out = 500)

png(file.path(output_folder, "conjugate_distribution_example.png"),
    width = 1400, height = 650, res = 150)
par(mfrow = c(1, 2))

# The observed sample's mean/SD differ from mapped prediction uncertainty.
hist(field_stock, col = "grey80", main = "Local field stocks",
     xlab = "Carbon stock (kg C/m2)")
abline(v = mean(field_stock), col = "orange", lwd = 2)
legend("topright", bty = "n", cex = 0.8,
       legend = c(paste("Mean:", round(mean(field_stock), 2)),
                  paste("SD:", round(sd(field_stock), 2))))

plot(x, dnorm(x, mu[1], sigma[1]), type = "l", col = "blue", lwd = 2,
     ylim = c(0, max(dnorm(mu, mu, sigma))),
     main = "Prior, local estimate and posterior", xlab = "Log carbon stock", ylab = "Density")
lines(x, dnorm(x, mu[2], sigma[2]), col = "orange", lwd = 2)
lines(x, dnorm(x, mu[3], sigma[3]), col = "darkgreen", lwd = 2)
legend("topright", legend = c("Prior", "Local likelihood (normalized)", "Posterior"),
       col = c("blue", "orange", "darkgreen"), lty = 1, lwd = 2, bty = "n", cex = 0.8)
dev.off()

# ----------------
# Step 4 - Plot and save the updated maps
# ----------------

plot(post_stock, main = "Area-based posterior mean (kg C/m2)")
plot(post_stock_idw, main = "IDW posterior mean (kg C/m2)")

# SD is pointwise uncertainty; it does not give uncertainty in an area total.
names(post_stock) <- "posterior_mean"
names(post_stock_sd) <- "posterior_sd"
names(post_stock_idw) <- "posterior_mean"
names(post_stock_sd_idw) <- "posterior_sd"
dir.create(file.path(output_folder, "area_based"), showWarnings = FALSE)
dir.create(file.path(output_folder, "idw"), showWarnings = FALSE)
writeRaster(post_stock, file.path(output_folder, "area_based/posterior_mean.tif"), overwrite = TRUE)
writeRaster(post_stock_sd, file.path(output_folder, "area_based/posterior_sd.tif"), overwrite = TRUE)
writeRaster(post_stock_idw, file.path(output_folder, "idw/posterior_mean.tif"), overwrite = TRUE)
writeRaster(post_stock_sd_idw, file.path(output_folder, "idw/posterior_sd.tif"), overwrite = TRUE)



# ----------------
# Step 5 - Spatial Bayesian update with inlabru
# ----------------

library(INLA)
library(inlabru)
library(fmesher)
library(sf)

# Simple assumptions for this first test.
gp_range_m <- 20000    # Fixed Matérn range: 20 km, not estimated.
core_log_sd <- 0.001   # Nearly exact observations; numerical approximation to zero.

# 1. Attach the McMaster prior to the core locations.
gp_cores <- st_as_sf(core_points)
p <- terra::extract(prior, core_points, ID = FALSE)

gp_cores$log_stock <- log(field_stock)
gp_cores$prior_m <- p$log_mean
gp_cores$prior_s <- sqrt(p$log_variance)

# 2. Attach the same prior to the prediction points.
gp_grid <- st_as_sf(predicted_points)
p <- terra::extract(prior, predicted_points, ID = FALSE)

gp_grid$prior_m <- p$log_mean
gp_grid$prior_s <- sqrt(p$log_variance)
gp_grid <- gp_grid[
  is.finite(gp_grid$prior_m) & is.finite(gp_grid$prior_s), ]

# 3. Make a triangular mesh extending beyond the study boundary.
gp_boundary <- st_buffer(st_as_sf(shared$aoi), dist = gp_range_m)

gp_mesh <- fm_mesh_2d(
  loc = gp_cores,
  boundary = gp_boundary,
  max.edge = c(gp_range_m / 4, gp_range_m),
  cutoff = 1000
)

# 4. Define a spatial field with fixed range and unit SD.
# Multiplying this field by prior_s supplies McMaster's local uncertainty.
gp_spde <- inla.spde2.pcmatern(
  gp_mesh,
  prior.range = c(gp_range_m, NA),
  prior.sigma = c(1, NA)
)

# Model:
# log carbon = McMaster log mean + McMaster log SD * spatial field.
# The field has prior mean zero, so predictions start centred on McMaster.
gp_fit <- bru(
  ~ -1 + field(geometry, model = gp_spde),
  bru_obs(
    log_stock ~ prior_m + prior_s * field,
    family = "gaussian",
    data = gp_cores,
    control.family = list(hyper = list(prec = list(
      initial = log(1 / core_log_sd^2),
      fixed = TRUE
    )))
  ),
  options = list(num.threads = "1:1")
)

# 5. Predict carbon stock by transforming posterior draws.
# This gives the carbon-scale mean and SD directly.
gp_prediction <- predict(
  gp_fit,
  gp_grid,
  ~ exp(prior_m + prior_s * field),
  n.samples = 1000,
  seed = 1
)

# 6. Put the results onto our existing raster grid.
gp_grid$mean <- gp_prediction$mean
gp_grid$sd <- gp_prediction$sd

gp_mean <- terra::rasterize(
  vect(gp_grid), analysis_mask, field = "mean"
)
gp_sd <- terra::rasterize(
  vect(gp_grid), analysis_mask, field = "sd"
)

plot(
  c(prior_mean, gp_mean, gp_sd),
  main = c(
    "McMaster mean",
    "Spatial posterior mean",
    "Spatial posterior SD"
  )
)

