# Run after scripts/03_area.R, 04_spatial.R and 05_bayes.R.
source("R/spatial.R")

# Equal precision gives an equal-weight mean and halves the variance.
z <- normal_update(10, 2, 20, 2)
stopifnot(abs(z$mean - 15) < 1e-12, abs(z$variance - 2) < 1e-12)
# A very weak new estimate leaves the prior essentially unchanged; reverse too.
stopifnot(abs(normal_update(10, 1, 30, 1e6)$mean - 10) < 1e-8)
stopifnot(abs(normal_update(10, 1e6, 30, 1)$mean - 30) < 1e-8)
z <- normal_update(c(1, 3), 2, c(5, 7), 3)
stopifnot(nrow(z) == 2, all(z$sd < 2), all(z$mean > c(1, 3)), all(z$mean < c(5, 7)))
stopifnot(inherits(try(normal_update(0, 0, 1, 1), silent = TRUE), "try-error"))

# Equal distances average the inputs; a sampled point reproduces its value.
stopifnot(abs(idw_values(matrix(c(2, 2), nrow = 1), c(10, 20)) - 15) < 1e-12)
stopifnot(abs(idw_values(matrix(c(0, 2), nrow = 1), c(10, 20)) - 10) < 1e-12)

area <- readRDS("outputs/part2/area/area_analysis.rds")
actual_area <- as.numeric(sf::st_area(area$aoi))
stopifnot(abs(sum(area$cells$area_m2)/actual_area - 1) < 1e-7)
stopifnot(abs(sum(area$area_totals$area_km2)*1e6/actual_area - 1) < 1e-7)
stopifnot(all(lengths(sf::st_intersects(area$points, area$aoi)) > 0))
expected <- sum(area$cells$area_m2*area$cells$stock_0_30_kg_m2)/1000
stopifnot(abs(sum(area$area_totals$total_stock_0_30_t_C) - expected) < 1e-6)

spatial <- readRDS("outputs/part2/spatial/spatial_analysis.rds")
train <- spatial$simulation[spatial$simulation$split == "Training", ]
test <- spatial$simulation[spatial$simulation$split == "Spatial holdout", ]
stopifnot(!length(intersect(train$cell_id, test$cell_id)), nrow(train) > nrow(test), nrow(test) > 10)
stopifnot(max(train$x_m) < min(test$x_m))
stopifnot(all(grepl("SIMULATED", spatial$simulation$record_origin)))
stopifnot(all(is.finite(spatial$cells$model_stock_0_30_kg_m2)))
stopifnot(all(spatial$cells$idw_stock_0_30_kg_m2 >= min(spatial$cores$stock_0_30_kg_m2) - 1e-10))
stopifnot(all(spatial$cells$idw_stock_0_30_kg_m2 <= max(spatial$cores$stock_0_30_kg_m2) + 1e-10))

posterior <- readRDS("outputs/part3/bayesian_analysis.rds")$cells
stopifnot(all(posterior$posterior_sd_kg_m2 < posterior$prior_sd_kg_m2))
stopifnot(all(posterior$posterior_sd_kg_m2 < posterior$new_sd_kg_m2))
stopifnot(all(posterior$posterior_mean_kg_m2 >= pmin(posterior$prior_mean_kg_m2, posterior$new_mean_kg_m2)))
stopifnot(all(posterior$posterior_mean_kg_m2 <= pmax(posterior$prior_mean_kg_m2, posterior$new_mean_kg_m2)))
message("Spatial area, IDW, simulation split and Bayesian checks passed.")
