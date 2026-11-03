# Verify branch completeness, data lineage, zero handling and independent update algebra.
library(terra)
source("config.R")
source("helpers/workshop_gaussian.R")
suite_checks <- list()
suite_check <- function(label, value) {
  suite_checks[[label]] <<- data.frame(check = label, passed = isTRUE(value))
  if (!isTRUE(value)) stop("Full-suite verification failed: ", label)
}
field_suite <- read.csv("outputs/script3d_li_suite/field_observations.csv")
eligible_suite <- field_suite[field_suite$canfit,]
suite_check("All authorized Li observations retain distinct evidence labels",
  nrow(eligible_suite) == 8 && sum(eligible_suite$confirmed) == 2 && sum(!eligible_suite$confirmed) == 6 &&
    sum(eligible_suite$depth == 0) == 4 && all(eligible_suite$carbon[eligible_suite$depth == 0] == 0))
layers_suite <- read.csv("data/field/core_layers.csv")
expected_stock_suite <- vapply(seq_len(nrow(eligible_suite)), function(i) {
  d <- layers_suite[layers_suite$plot_id == eligible_suite$core_id[i],]
  sum(d$soc*d$bulk_density*(1-d$coarse_frag)*
    pmax(0,pmin(d$depth_to,eligible_suite$depth[i])-d$depth_from)/100)
}, numeric(1))
suite_check("Organic stocks independently close to raw layer integration",
  max(abs(expected_stock_suite-eligible_suite$carbon)) < 1e-10)
suite_methods <- c("raw_cores","area_based","idw","regression")
manifest_suite <- read.csv("outputs/script5/li_suite/completion_manifest.csv")
suite_check("All eight Li update branches exist", nrow(manifest_suite) == 8 &&
  all(manifest_suite$completed) && all(file.exists(manifest_suite$update)) &&
  setequal(paste(manifest_suite$variable,manifest_suite$method),
    as.vector(outer(c("depth","carbon"),suite_methods,paste))))

for (variable_suite in c("depth","carbon")) {
  base_suite <- file.path("outputs/script5/li_suite",variable_suite)
  reference_cells_suite <- NULL
  for (method_suite in suite_methods) {
    d <- read.csv(file.path(base_suite,method_suite,"predictions.csv"))
    if (is.null(reference_cells_suite)) reference_cells_suite <- d$cell
    suite_check(paste(variable_suite,method_suite,"uses common prior cells"), identical(d$cell,reference_cells_suite))
    suite_check(paste(variable_suite,method_suite,"finite means, valid conditional variance and bounds"),
      all(is.finite(d$posterior_mean)) && all(is.finite(d$posterior_sd)) && all(d$posterior_sd >= 0) &&
      all(d$posterior_sd^2 <= d$prior_sd^2+1e-7) &&
      all(d$lower95 <= d$posterior_mean & d$posterior_mean <= d$upper95))
    suite_check(paste(variable_suite,method_suite,"negative posterior flags match values"),
      all(d$negative_mean == as.integer(d$posterior_mean < 0)))
    s <- read.csv(file.path(base_suite,method_suite,"summary.csv"))
    suite_check(paste(variable_suite,method_suite,"all eight cores and four zero responses retained"),
      s$n_cores == 8 && s$zero_observations == 4)
    if (method_suite != "raw_cores") {
      applied <- d$update_applied == 1
      # Independent precision form, not the weighting implementation.
      variance <- 1/(1/d$prior_sd[applied]^2+1/s$working_error_sd^2)
      mu <- variance*(d$prior_mean[applied]/d$prior_sd[applied]^2+
        d$local_estimate[applied]/s$working_error_sd^2)
      suite_check(paste(variable_suite,method_suite,"independent precision-form closure"),
        max(abs(mu-d$posterior_mean[applied])) < 1e-8 &&
          max(abs(variance-d$posterior_sd[applied]^2)) < 1e-7)
      retained <- !applied
      suite_check(paste(variable_suite,method_suite,"unsupported local evidence retains prior"),
        all(d$posterior_mean[retained] == d$prior_mean[retained]) &&
          all(d$data_weight[retained] == 0) && all(d$local_negative_or_missing[retained] == 1))
      suite_check(paste(variable_suite,method_suite,"no negative local prediction used as evidence"),
        all(d$local_estimate[applied] >= 0))
    }
    curve_suite <- rast(file.path("outputs/script5/li_suite/depth",method_suite,"curve_derived_carbon.tif"))
    depth_suite <- rast(file.path("outputs/script5/li_suite/depth",method_suite,"posterior_maps.tif"))$posterior_mean
    finite_curve <- is.finite(values(curve_suite[[1]],mat=FALSE))
    dd <- values(depth_suite,mat=FALSE)[finite_curve]
    suite_check(paste(method_suite,"curve mapping stays within measured depth support",variable_suite),
      all(dd >= 0 & dd <= max(eligible_suite$depth)))
  }
  cv_suite <- read.csv(file.path(base_suite,"nested_leave_one_core_out.csv"))
  suite_check(paste(variable_suite,"all four methods have eight complete outer holdouts"),
    nrow(cv_suite) == 32 && all(table(cv_suite$method) == 8) &&
      !anyDuplicated(cv_suite[c("method","core_id")]) && all(is.finite(cv_suite$predicted)))
  fold_suite <- read.csv(file.path(base_suite,"nested_fold_choices.csv"))
  suite_check(paste(variable_suite,"nested model choices recorded per held-out core"),
    nrow(fold_suite) == 8 && setequal(fold_suite$core_id,eligible_suite$core_id))
  prior_audit_suite <- read.csv(file.path(base_suite,"core_prior_extraction_audit.csv"))
  suite_check(paste(variable_suite,"paired prior fallback bounded and all core priors finite"),
    all(prior_audit_suite$source_distance_m <= 500) && all(is.finite(prior_audit_suite$used_mean)))
}
exact_zero_suite <- workshop_gaussian_gp(matrix(c(0,0),1),0,20,5,matrix(c(0,0),1),20,5)
suite_check("Original-scale GP exactly conditions on true zero without an offset",
  abs(exact_zero_suite$posterior_mean) < 1e-10 && exact_zero_suite$posterior_sd < 1e-8)
zero_prior_suite <- workshop_gaussian_gp(matrix(c(0,0),1),0,0,5,matrix(c(1,0),1),0,5)
suite_check("Valid zero-valued prior means remain usable", is.finite(zero_prior_suite$posterior_mean))
curve_summary_suite <- read.csv("outputs/script3d_li_suite/carbon_curve/model_summary.csv")
suite_check("Zero-depth profiles do not inflate carbon-curve residual degrees of freedom",
  curve_summary_suite$n_positive_profiles_fitted == 4 && curve_summary_suite$residual_df == 3)
for (method_suite in suite_methods) {
  path_suite <- if (method_suite == "raw_cores") "raw_core_maps.tif" else "posterior_maps.tif"
  suite_check(paste("Sothe",method_suite,"update exists"),file.exists(file.path(
    "outputs/script5/sothe/stock_0_30",method_suite,path_suite)))
}
regression_checks_suite <- read.csv("outputs/script5/sothe/stock_0_30/regression/verification.csv")
suite_check("Fourth Sothe regression update numerical checks pass",all(regression_checks_suite$passed))
combined_manifest_suite <- rbind(
  data.frame(variable = paste0("li_",manifest_suite$variable),method = manifest_suite$method,
    file = manifest_suite$update, completed = manifest_suite$completed),
  do.call(rbind,lapply(suite_methods,function(m) data.frame(variable = "sothe_s30",method = m,
    file = file.path("outputs/script5/sothe/stock_0_30",m,
      if(m == "raw_cores") "raw_core_maps.tif" else "posterior_maps.tif"),completed = TRUE))))
suite_check("Full comparison has twelve completed update outputs",
  nrow(combined_manifest_suite) == 12 && all(file.exists(combined_manifest_suite$file)) &&
    all(table(combined_manifest_suite$variable) == 4))
full_suite_checks <- do.call(rbind,suite_checks)
dir.create("outputs/verification",recursive=TRUE,showWarnings=FALSE)
write.csv(full_suite_checks,"outputs/verification/full_suite_checks.csv",row.names=FALSE)
write.csv(combined_manifest_suite,"outputs/verification/completed_updates.csv",row.names=FALSE)
message("Full-suite verification passed: ",nrow(full_suite_checks)," checks.")
