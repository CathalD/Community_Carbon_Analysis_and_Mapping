#!/usr/bin/env Rscript

# Outcome-blind external-neighbour rankings; evaluate matched 0-30 cm stock.
if (!file.exists("config/route_l.R")) stop("Run from the project root.")
source("R/covariates/alphaearth_diagnostics.R")

external <- utils::read.csv(
  "data/derived/covariate_analysis/external_core_covariates_riley_canpeat.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
community <- utils::read.csv(
  "data/derived/covariate_analysis/community_core_covariates.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
folds <- utils::read.csv(
  "outputs/covariate_analysis/part5_external_model/external_resampling_assignments.csv",
  stringsAsFactors = FALSE
)
cfg <- dget("config/covariate_analysis.R")
env_vars <- cfg$local_similarity_experiment$interpretable_environment$variables
ae_vars <- alphaearth_columns(external)
h_ae <- cfg$local_similarity_experiment$alphaearth$primary_bandwidth_observed
stock_col <- "carbon_stock_0_30cm_kg_m2"

external$outer_fold <- folds$outer_fold[match(external$core_uid, folds$core_uid)]
if (any(is.na(external$outer_fold[external$carbon_complete_0_30cm]))) {
  stop("Missing locked spatial/study fold for an external complete stock.")
}
external_ae <- l2_normalize_rows(external, ae_vars)
external_keep <- external$carbon_complete_0_30cm &
  is.finite(external[[stock_col]]) &
  stats::complete.cases(external[, env_vars, drop = FALSE]) &
  external_ae$valid
external_eligible <- external[external_keep, , drop = FALSE]
external_eligible$site_key <- paste(
  external_eligible$database_code,
  sprintf("%.6f", external_eligible$longitude),
  sprintf("%.6f", external_eligible$latitude), sep = "|"
)
site_groups <- split(seq_len(nrow(external_eligible)), external_eligible$site_key)
if (any(vapply(site_groups, function(i) length(unique(external_eligible$outer_fold[i])), integer(1)) != 1L)) {
  stop("Profiles from one external site cross held-out folds.")
}
sites <- external_eligible[vapply(site_groups, `[`, integer(1), 1L), , drop = FALSE]
sites$n_profiles_at_site <- lengths(site_groups)
sites$response_kg_m2 <- vapply(
  site_groups, function(i) mean(external_eligible[[stock_col]][i]), numeric(1)
)
sites$target_id <- sites$site_key
rownames(sites) <- NULL
if (anyDuplicated(sites$site_key)) stop("Site aggregation failed.")
coordinate_key <- paste(sprintf("%.6f", sites$longitude),
                        sprintf("%.6f", sites$latitude), sep = "|")
if (anyDuplicated(coordinate_key)) {
  stop("External sources share rounded coordinates; cross-source test would need site grouping.")
}

community_ae <- l2_normalize_rows(community, ae_vars)
community_keep <- community$coverage_0_30_fraction >= 1 - 1e-8 &
  is.finite(community$measured_stock_0_30_kg_m2) &
  stats::complete.cases(community[, env_vars, drop = FALSE]) &
  community_ae$valid
community_test <- community[community_keep, , drop = FALSE]
community_test$response_kg_m2 <- community_test$measured_stock_0_30_kg_m2
community_test$target_id <- community_test$plot_id
if (nrow(community_test) != 6L) stop("Expected six complete community 0-30 cm profiles.")

as_numeric_matrix <- function(data, columns) {
  z <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(z) <- "double"
  z
}
haversine_matrix_km <- function(target, donor) {
  rad <- pi / 180
  dlat <- outer(target$latitude, donor$latitude, "-") * rad
  dlon <- outer(target$longitude, donor$longitude, "-") * rad
  a <- sin(dlat / 2)^2 +
    outer(cos(target$latitude * rad), cos(donor$latitude * rad), "*") * sin(dlon / 2)^2
  a[a < 0] <- 0
  a[a > 1] <- 1
  2 * 6371.0088 * asin(sqrt(a))
}
score_matrices <- function(target, donor) {
  donor_env <- as_numeric_matrix(donor, env_vars)
  target_env <- as_numeric_matrix(target, env_vars)
  centre <- apply(donor_env, 2, stats::median)
  scale <- apply(donor_env, 2, stats::mad)
  fallback <- apply(donor_env, 2, stats::sd)
  use_fallback <- !is.finite(scale) | scale <= 0
  scale[use_fallback] <- fallback[use_fallback]
  scale[!is.finite(scale) | scale <= 0] <- 1
  donor_z <- sweep(sweep(donor_env, 2, centre, "-"), 2, scale, "/")
  target_z <- sweep(sweep(target_env, 2, centre, "-"), 2, scale, "/")
  env <- matrix(0, nrow(target), nrow(donor))
  for (j in seq_along(env_vars)) {
    env <- env + exp(-0.5 * outer(target_z[, j], donor_z[, j], "-")^2) / length(env_vars)
  }
  ae_t <- l2_normalize_rows(target, ae_vars)
  ae_d <- l2_normalize_rows(donor, ae_vars)
  if (!all(ae_t$valid) || !all(ae_d$valid)) stop("Invalid AlphaEarth vector in comparison pool.")
  cosine <- tcrossprod(ae_t$normalized, ae_d$normalized)
  cosine[cosine < -1] <- -1
  cosine[cosine > 1] <- 1
  squared_chord <- 2 - 2 * cosine
  squared_chord[squared_chord < 0] <- 0
  chord <- sqrt(squared_chord)
  ae <- exp(-0.5 * (chord / h_ae)^2)
  list(environment = env, alphaearth = ae,
       geographic_km = haversine_matrix_km(target, donor))
}

all_predictions <- list()
all_associations <- list()
community_neighbours <- list()
prediction_index <- 0L
association_index <- 0L
neighbour_index <- 0L
evaluate_pool <- function(target, donor, cohort, source_code, assessment_fold = NA_integer_) {
  scores <- score_matrices(target, donor)
  donor_stock <- donor$response_kg_m2
  n_donor <- nrow(donor)
  if (n_donor < 50L) stop("Fewer than 50 available donor sites.")
  for (i in seq_len(nrow(target))) {
    observed <- target$response_kg_m2[i]
    donor_abs_error <- abs(donor_stock - observed)
    env_rank <- rank(-scores$environment[i, ], ties.method = "average")
    ae_rank <- rank(-scores$alphaearth[i, ], ties.method = "average")
    combined_rank <- env_rank + ae_rank
    orderings <- list(
      geographic = order(scores$geographic_km[i, ], donor$site_key),
      environment = order(-scores$environment[i, ], donor$site_key),
      alphaearth = order(-scores$alphaearth[i, ], donor$site_key),
      combined_rank = order(combined_rank, donor$site_key)
    )
    associations <- c(
      geographic = suppressWarnings(stats::cor(-scores$geographic_km[i, ],
                                                -donor_abs_error, method = "spearman")),
      environment = suppressWarnings(stats::cor(scores$environment[i, ],
                                                 -donor_abs_error, method = "spearman")),
      alphaearth = suppressWarnings(stats::cor(scores$alphaearth[i, ],
                                                -donor_abs_error, method = "spearman"))
    )
    for (method in names(associations)) {
      association_index <<- association_index + 1L
      all_associations[[association_index]] <<- data.frame(
        cohort = cohort, source = source_code, target_id = target$target_id[i],
        assessment_fold = assessment_fold, method = method,
        spearman_similarity_vs_negative_donor_error = unname(associations[method]),
        n_donor_sites = n_donor, stringsAsFactors = FALSE
      )
    }
    prediction_index <<- prediction_index + 1L
    all_predictions[[prediction_index]] <<- data.frame(
      cohort = cohort, source = source_code, target_id = target$target_id[i],
      assessment_fold = assessment_fold, method = "source_mean", k = NA_integer_,
      observed_kg_m2 = observed, predicted_kg_m2 = mean(donor_stock),
      mean_absolute_donor_difference_kg_m2 = mean(donor_abs_error),
      all_donor_mean_absolute_difference_kg_m2 = mean(donor_abs_error),
      median_donor_distance_km = stats::median(scores$geographic_km[i, ]),
      n_donor_sites = n_donor, stringsAsFactors = FALSE
    )
    for (k in c(25L, 50L)) for (method in names(orderings)) {
      chosen <- orderings[[method]][seq_len(k)]
      prediction_index <<- prediction_index + 1L
      all_predictions[[prediction_index]] <<- data.frame(
        cohort = cohort, source = source_code, target_id = target$target_id[i],
        assessment_fold = assessment_fold, method = method, k = k,
        observed_kg_m2 = observed, predicted_kg_m2 = mean(donor_stock[chosen]),
        mean_absolute_donor_difference_kg_m2 = mean(donor_abs_error[chosen]),
        all_donor_mean_absolute_difference_kg_m2 = mean(donor_abs_error),
        median_donor_distance_km = stats::median(scores$geographic_km[i, chosen]),
        n_donor_sites = n_donor, stringsAsFactors = FALSE
      )
      if (cohort == "community" && k == 25L && method != "geographic") {
        top <- chosen
        neighbour_index <<- neighbour_index + 1L
        community_neighbours[[neighbour_index]] <<- data.frame(
          target_id = target$target_id[i], source = source_code, method = method,
          rank = seq_len(25L), site_key = donor$site_key[top],
          representative_core_uid = donor$core_uid[top],
          study_id = donor$study_id[top],
          n_profiles_at_site = donor$n_profiles_at_site[top],
          ecosystem = donor$ecosystem[top],
          latitude = donor$latitude[top], longitude = donor$longitude[top],
          distance_km = scores$geographic_km[i, top],
          environmental_similarity = scores$environment[i, top],
          alphaearth_similarity = scores$alphaearth[i, top],
          observed_community_kg_m2 = observed,
          external_site_mean_kg_m2 = donor_stock[top],
          absolute_stock_difference_kg_m2 = donor_abs_error[top],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

for (source_code in c("CANPEAT", "RILEY")) {
  source_sites <- sites[sites$database_code == source_code, , drop = FALSE]
  evaluate_pool(community_test, source_sites, "community", source_code)
  for (fold in sort(unique(source_sites$outer_fold))) {
    target <- source_sites[source_sites$outer_fold == fold, , drop = FALSE]
    donor <- source_sites[source_sites$outer_fold != fold, , drop = FALSE]
    evaluate_pool(target, donor, "external_spatial_study_holdout", source_code, fold)
  }
}

# A harder transfer check: an entire source is unseen by the donor pool.
# This does not substitute for the community evaluation, but exposes source shift.
evaluate_pool(
  sites[sites$database_code == "RILEY", , drop = FALSE],
  sites[sites$database_code == "CANPEAT", , drop = FALSE],
  "cross_source_external", "CANPEAT_to_RILEY"
)
evaluate_pool(
  sites[sites$database_code == "CANPEAT", , drop = FALSE],
  sites[sites$database_code == "RILEY", , drop = FALSE],
  "cross_source_external", "RILEY_to_CANPEAT"
)

predictions <- do.call(rbind, all_predictions)
predictions$error_kg_m2 <- predictions$predicted_kg_m2 - predictions$observed_kg_m2
predictions$absolute_error_kg_m2 <- abs(predictions$error_kg_m2)
associations <- do.call(rbind, all_associations)
neighbours <- do.call(rbind, community_neighbours)
rownames(predictions) <- rownames(associations) <- rownames(neighbours) <- NULL

summarize <- function(data) {
  groups <- split(data, interaction(data$cohort, data$source, data$method,
                                    ifelse(is.na(data$k), "all", data$k), drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(z) data.frame(
    cohort = z$cohort[1], source = z$source[1], method = z$method[1], k = z$k[1],
    n_targets = nrow(z), MAE_kg_m2 = mean(z$absolute_error_kg_m2),
    RMSE_kg_m2 = sqrt(mean(z$error_kg_m2^2)), bias_kg_m2 = mean(z$error_kg_m2),
    mean_absolute_donor_difference_kg_m2 = mean(z$mean_absolute_donor_difference_kg_m2),
    mean_all_donor_difference_kg_m2 = mean(z$all_donor_mean_absolute_difference_kg_m2),
    median_selected_donor_distance_km = stats::median(z$median_donor_distance_km),
    stringsAsFactors = FALSE
  )))
  rownames(out) <- NULL
  out[order(out$cohort, out$source, out$method, out$k), ]
}
summary <- summarize(predictions)

# Equal source weighting prevents CanPEAT's larger site count from silently
# overwhelming Riley in a pooled community prediction.
community_predictions <- predictions[predictions$cohort == "community", ]
balance_key <- paste(community_predictions$target_id, community_predictions$method,
                     community_predictions$k, sep = "|")
balanced <- do.call(rbind, lapply(split(seq_len(nrow(community_predictions)), balance_key), function(i) {
  z <- community_predictions[i, ]
  if (nrow(z) != 2L || !setequal(z$source, c("CANPEAT", "RILEY"))) {
    stop("Source-balanced row lacks one of the two external databases.")
  }
  y <- z[1, , drop = FALSE]
  y$source <- "EQUAL_SOURCE_BLEND"
  y$predicted_kg_m2 <- mean(z$predicted_kg_m2)
  y$error_kg_m2 <- y$predicted_kg_m2 - y$observed_kg_m2
  y$absolute_error_kg_m2 <- abs(y$error_kg_m2)
  y$mean_absolute_donor_difference_kg_m2 <- mean(z$mean_absolute_donor_difference_kg_m2)
  y$all_donor_mean_absolute_difference_kg_m2 <- mean(z$all_donor_mean_absolute_difference_kg_m2)
  y$median_donor_distance_km <- stats::median(z$median_donor_distance_km)
  y$n_donor_sites <- sum(z$n_donor_sites)
  y
}))
rownames(balanced) <- NULL
summary <- rbind(summary, summarize(balanced))

association_summary <- do.call(rbind, lapply(
  split(associations, interaction(associations$cohort, associations$source,
                                  associations$method, drop = TRUE)),
  function(z) data.frame(
    cohort = z$cohort[1], source = z$source[1], method = z$method[1],
    n_targets = nrow(z),
    median_spearman = stats::median(z$spearman_similarity_vs_negative_donor_error,
                                    na.rm = TRUE),
    positive_fraction = mean(z$spearman_similarity_vs_negative_donor_error > 0,
                             na.rm = TRUE), stringsAsFactors = FALSE
  )
))
rownames(association_summary) <- NULL

out_dir <- file.path("outputs", "covariate_analysis", "part15_external_similarity_transfer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sites[, c("site_key", "core_uid", "database_code", "study_id",
                           "latitude", "longitude", "ecosystem", "n_profiles_at_site",
                           "outer_fold", "response_kg_m2")],
                 file.path(out_dir, "external_unique_sites.csv"), row.names = FALSE)
utils::write.csv(predictions, file.path(out_dir, "source_specific_predictions.csv"), row.names = FALSE)
utils::write.csv(balanced, file.path(out_dir, "community_equal_source_predictions.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "prediction_summary.csv"), row.names = FALSE)
utils::write.csv(associations, file.path(out_dir, "similarity_stock_associations.csv"), row.names = FALSE)
utils::write.csv(association_summary,
                 file.path(out_dir, "similarity_stock_association_summary.csv"), row.names = FALSE)
utils::write.csv(neighbours[neighbours$rank <= 10L, ],
                 file.path(out_dir, "community_top10_external_sites.csv"), row.names = FALSE)
utils::write.csv(neighbours, file.path(out_dir, "community_top25_external_sites.csv"), row.names = FALSE)
writeLines(c(
  "# Riley and CanPEAT similarity transfer protocol",
  "",
  "Only complete, measured 0-30 cm carbon stocks enter the comparison. External profiles with missing five-variable or AlphaEarth covariates are excluded from every candidate so results share one analysis pool.",
  "External profiles at the same source/coordinate are averaged to one site before neighbour selection. All profiles from a site share one locked spatial/study fold.",
  "Environmental similarity is the equal mean of five Gaussian similarities on donor-site median/MAD scales. AlphaEarth uses L2-normalized 64-band vectors and chord distance. The combined rule adds their within-pool ranks equally. No local carbon outcome selects features, bandwidth, neighbour count or source mix.",
  "Predictions are unweighted means of the 25 or 50 most similar donor sites; controls are all-source site mean and the 25/50 geographically nearest external sites. Community comparisons are shown separately for Riley and CanPEAT and in an equal-source blend.",
  "Within-source external checks use the existing connected spatial/study outer folds, with every target site absent from its donor pool. Cross-source checks train on all sites of one database and test on every site of the other; no exact source-coordinate overlap was detected. The community test contains only six complete profiles, so differences are descriptive and not an independent model-selection result.",
  "Donor-error comparisons and Spearman associations ask whether more similar sites have more similar measured stock. No p-values are reported because site relationships, folds and six community targets are not independent replicates.",
  "These are fixed 0-30 cm mineral-plus-organic profile stocks. They do not validate the Route L organic contact stock or full variable-depth total. No raster or authoritative map is changed."
), file.path(out_dir, "METHODS.md"))
cat("Eligible external profiles:", nrow(external_eligible),
    "; unique sites:", nrow(sites), "; community targets:", nrow(community_test), "\n")
print(summary[summary$cohort == "community" &
                (summary$method == "source_mean" | summary$k == 25L),
              c("source", "method", "k", "n_targets", "MAE_kg_m2",
                "mean_absolute_donor_difference_kg_m2")], row.names = FALSE)
cat("Wrote:", out_dir, "\n")
