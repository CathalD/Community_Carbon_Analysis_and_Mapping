# Inputs: site table and eligible original-study targets from historical script 79.
# New exploratory log-carbon donor comparison; no GP and no posterior uncertainty.
stopifnot(is.finite(environment_lengthscale), environment_lengthscale > 0,
          top_donors_to_export >= 1L)
log_rows <- list()
top_rows <- list()
scaling_rows <- list()
log_index <- top_index <- scaling_index <- 0L

log_donor_pool <- function(target, donor, cohort, source, fold = NA_integer_) {
  stopifnot(nrow(donor) >= 2L, !anyDuplicated(donor$site_key),
            !any(target$site_key %in% donor$site_key),
            all(donor$response_kg_m2 > 0), all(target$response_kg_m2 > 0))
  x <- as_numeric_matrix(donor, env_vars)
  xt <- as_numeric_matrix(target, env_vars)
  centre <- apply(x, 2, median)
  spread <- apply(x, 2, mad)
  fallback <- apply(x, 2, sd)
  bad <- !is.finite(spread) | spread <= 0
  spread[bad] <- fallback[bad]
  spread[!is.finite(spread) | spread <= 0] <- 1
  z <- sweep(sweep(x, 2, centre, "-"), 2, spread, "/")
  zt <- sweep(sweep(xt, 2, centre, "-"), 2, spread, "/")
  d2 <- matrix(0, nrow(target), nrow(donor))
  for (j in seq_along(env_vars)) d2 <- d2 + outer(zt[, j], z[, j], "-")^2
  env_log_s <- -0.5 * d2 / environment_lengthscale^2
  ae_d <- l2_normalize_rows(donor, ae_vars)$normalized
  ae_t <- l2_normalize_rows(target, ae_vars)$normalized
  # Estimate bandwidth only from donors belonging to this training pool.
  train_d2 <- pmax(2 - 2 * tcrossprod(ae_d), 0)
  diag(train_d2) <- Inf
  h <- median(sqrt(apply(train_d2, 1, min)))
  if (!is.finite(h) || h <= 0) stop("Degenerate training-pool AlphaEarth bandwidth.")
  ae_d2 <- pmax(2 - 2 * tcrossprod(ae_t, ae_d), 0)
  ae_log_s <- -0.5 * ae_d2 / h^2
  scaling_index <<- scaling_index + 1L
  scaling_rows[[scaling_index]] <<- data.frame(
    cohort = cohort, source = source, fold = fold, predictor = env_vars,
    training_median = centre, training_scale = spread,
    environment_lengthscale = environment_lengthscale,
    training_alphaearth_bandwidth = h, training_sites = nrow(donor))
  y <- log(donor$response_kg_m2)
  for (i in seq_len(nrow(target))) {
    scores <- list(source_log_mean = rep(0, nrow(donor)),
                   environmental_kernel = env_log_s[i, ],
                   alphaearth_kernel = ae_log_s[i, ])
    for (method in names(scores)) {
      ls <- scores[[method]]
      # Stable normalized similarities even when all absolute weights are tiny.
      w <- exp(ls - max(ls)); w <- w / sum(w)
      estimate <- sum(w * y)
      donor_spread <- sqrt(sum(w * (y - estimate)^2))
      log_index <<- log_index + 1L
      log_rows[[log_index]] <<- data.frame(
        cohort = cohort, source = source, fold = fold,
        target_id = target$target_id[i], method = method,
        observed_kg_m2 = target$response_kg_m2[i], log_estimate = estimate,
        geometric_stock_estimate_kg_m2 = exp(estimate),
        donor_log_spread = donor_spread, effective_donor_sites = 1 / sum(w^2),
        maximum_similarity = exp(max(ls)), maximum_log_similarity = max(ls),
        nearest_environment_distance = sqrt(min(d2[i, ])),
        n_environment_predictors_outside_training_range = sum(
          xt[i, ] < apply(x, 2, min) | xt[i, ] > apply(x, 2, max)),
        donor_sites = nrow(donor), stringsAsFactors = FALSE)
      if (cohort == "community" && method != "source_log_mean") {
        chosen <- head(order(-w, donor$site_key), top_donors_to_export)
        top_index <<- top_index + 1L
        top_rows[[top_index]] <<- data.frame(
          target_id = target$target_id[i], source = source, method = method,
          rank = seq_along(chosen), core_uid = donor$core_uid[chosen],
          study_id = donor$study_id[chosen], ecosystem = donor$ecosystem[chosen],
          latitude = donor$latitude[chosen], longitude = donor$longitude[chosen],
          site_stock_kg_m2 = donor$response_kg_m2[chosen],
          normalized_donor_weight = w[chosen], similarity = exp(ls[chosen]),
          stringsAsFactors = FALSE)
      }
    }
  }
}

# Community responses remain absent from every fit and scale estimate.
community_log <- community_test
community_log$site_key <- paste0("COMMUNITY|", community_log$target_id)
for (source in c("CANPEAT", "RILEY")) {
  donor <- sites[sites$database_code == source, , drop = FALSE]
  log_donor_pool(community_log, donor, "community", source)
  for (fold in sort(unique(donor$outer_fold))) {
    log_donor_pool(donor[donor$outer_fold == fold, , drop = FALSE],
                   donor[donor$outer_fold != fold, , drop = FALSE],
                   "external_spatial_study_holdout", source, fold)
  }
}
for (source in c("CANPEAT", "RILEY")) {
  log_donor_pool(sites[sites$database_code != source, , drop = FALSE],
                 sites[sites$database_code == source, , drop = FALSE],
                 "cross_source_external", source)
}
log_predictions <- do.call(rbind, log_rows)
log_predictions$error_kg_m2 <- log_predictions$geometric_stock_estimate_kg_m2 -
  log_predictions$observed_kg_m2
keys <- interaction(log_predictions$cohort, log_predictions$source,
                    log_predictions$method, drop = TRUE)
log_metrics <- do.call(rbind, lapply(split(log_predictions, keys), function(d) {
  data.frame(cohort = d$cohort[1], source = d$source[1], method = d$method[1],
             n = nrow(d), MAE_kg_m2 = mean(abs(d$error_kg_m2)),
             RMSE_kg_m2 = sqrt(mean(d$error_kg_m2^2)),
             bias_kg_m2 = mean(d$error_kg_m2),
             RMSE_log = sqrt(mean((d$log_estimate - log(d$observed_kg_m2))^2)))
}))
log_out <- "outputs/part5_log_donors"
dir.create(log_out, recursive = TRUE, showWarnings = FALSE)
write.csv(log_predictions, file.path(log_out, "predictions.csv"), row.names = FALSE)
write.csv(log_metrics, file.path(log_out, "metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, top_rows), file.path(log_out, "community_top_donors.csv"), row.names = FALSE)
write.csv(do.call(rbind, scaling_rows), file.path(log_out, "training_pool_scaling.csv"), row.names = FALSE)
write.csv(sites, file.path(log_out, "eligible_external_sites.csv"), row.names = FALSE)
write.csv(community_test, file.path(log_out, "original_study_complete_targets.csv"), row.names = FALSE)
png(file.path(log_out, "community_comparison.png"), width = 1500, height = 900, res = 150)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (source in c("CANPEAT", "RILEY")) {
  for (method in c("source_log_mean", "environmental_kernel", "alphaearth_kernel")) {
    d <- log_predictions[log_predictions$cohort == "community" &
                  log_predictions$source == source & log_predictions$method == method, ]
    limits <- range(c(d$observed_kg_m2, d$geometric_stock_estimate_kg_m2))
    plot(d$observed_kg_m2, d$geometric_stock_estimate_kg_m2,
         xlim = limits, ylim = limits, pch = 19,
         xlab = "Observed stock (kg C/m2)", ylab = "Geometric stock estimate (kg C/m2)",
         main = paste(source, method, sep = "\n"))
    abline(0, 1, col = "grey50", lty = 2)
  }
}
dev.off()
print(log_metrics[log_metrics$cohort == "community", ], row.names = FALSE)
