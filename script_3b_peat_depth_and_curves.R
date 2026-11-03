# Additional workshop-style local models, before using Li priors.
# Sampled depth is not automatically peat depth. Only explicit field contacts
# establish a complete organic profile; other cores remain in Script 1 / Sothe.
library(readr)
library(dplyr)
library(ggplot2)
source("config.R")
out <- "outputs/script3b_li_local"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
layers <- read_csv("outputs/script1/field_layers.csv", show_col_types = FALSE)
cores <- read_csv("outputs/script1/mapping_data.csv", show_col_types = FALSE)
contacts <- read_csv("data/field/organic_mineral_contacts.csv", show_col_types = FALSE)
stopifnot(!anyDuplicated(contacts$plot_id), all(contacts$organic_horizon_depth_cm > 0))
observations <- lapply(seq_len(nrow(cores)), function(i) {
  id <- cores[["Core ID"]][i]; d <- layers[layers[["Core ID"]] == id, ]
  contact <- contacts[contacts$plot_id == id, ]
  valid <- nrow(contact) == 1 && contact$organic_horizon_depth_cm <= max(d[["Bottom depth (cm)"]]) + 1e-8
  depth <- if (valid) contact$organic_horizon_depth_cm else NA_real_
  stock <- if (valid) sum(d[["Stock, full thickness (kg C/m²)"]] /
    (d[["Bottom depth (cm)"]] - d[["Top depth (cm)"]]) *
    pmax(0, pmin(d[["Bottom depth (cm)"]], depth) - d[["Top depth (cm)"]])) else NA_real_
  data.frame(core_id = id, longitude = cores$Longitude[i], latitude = cores$Latitude[i],
    peat_depth_cm = depth, organic_carbon_kg_m2 = stock,
    sampled_bottom_cm = max(d[["Bottom depth (cm)"]]), eligible_for_li = valid,
    evidence = if (valid) contact$contact_source else "No explicit peat/mineral contact supplied")
})
all <- bind_rows(observations)
write_csv(all, file.path(out, "contact_eligibility.csv"))
local <- all[all$eligible_for_li, ]
if (nrow(local) < 2) stop("At least two confirmed peat contacts needed for exploratory local models")
write_csv(local, file.path(out, "li_field_observations.csv"))

# Script 3 mean-only baseline. Two contacts do not support covariate regression.
depth_model <- lm(peat_depth_cm ~ 1, local)
depth_cv <- vapply(seq_len(nrow(local)), function(i)
  as.numeric(predict(lm(peat_depth_cm ~ 1, local[-i,,drop=FALSE]), local[i,,drop=FALSE])), numeric(1))
write_csv(data.frame(core_id = local$core_id, observed_cm = local$peat_depth_cm,
  predicted_cm = depth_cv, error_cm = depth_cv-local$peat_depth_cm), file.path(out, "depth_leave_one_core_out.csv"))
saveRDS(depth_model, file.path(out, "depth_mean_model.rds"))

# Cumulative carbon with depth (not a rate through time). Fit one endpoint per
# independent core, through zero; do not count cumulative sections as independent.
carbon_curve <- lm(organic_carbon_kg_m2 ~ 0 + peat_depth_cm, local)
curve_slope <- unname(coef(carbon_curve)[1])
curve_cv <- vapply(seq_len(nrow(local)), function(i) {
  fit <- lm(organic_carbon_kg_m2 ~ 0 + peat_depth_cm, local[-i,,drop=FALSE])
  as.numeric(predict(fit, local[i,,drop=FALSE]))
}, numeric(1))
write_csv(data.frame(core_id = local$core_id, observed = local$organic_carbon_kg_m2,
  predicted = curve_cv, error = curve_cv-local$organic_carbon_kg_m2), file.path(out, "carbon_curve_leave_one_core_out.csv"))
saveRDS(carbon_curve, file.path(out, "carbon_depth_curve.rds"))
supported <- seq(0, max(local$peat_depth_cm), length.out = 101)
write_csv(data.frame(depth_cm = supported, carbon_kg_m2 = curve_slope*supported),
  file.path(out, "supported_cumulative_curve.csv"))
knots <- bind_rows(lapply(local$core_id, function(id) {
  d <- layers[layers[["Core ID"]] == id, ]
  limit <- local$peat_depth_cm[local$core_id == id]
  d <- d[d[["Top depth (cm)"]] < limit, ]
  stock <- d[["Stock, full thickness (kg C/m²)"]] /
    (d[["Bottom depth (cm)"]] - d[["Top depth (cm)"]]) *
    (pmin(d[["Bottom depth (cm)"]], limit)-d[["Top depth (cm)"]])
  data.frame(core_id = id, depth_cm = c(0, pmin(d[["Bottom depth (cm)"]], limit)),
    cumulative_stock_kg_m2 = c(0, cumsum(stock)))
}))
write_csv(knots, file.path(out, "measured_organic_cumulative_curves.csv"))
p <- ggplot(knots, aes(depth_cm, cumulative_stock_kg_m2, colour = core_id)) +
  geom_line() + geom_point() + geom_abline(intercept = 0, slope = curve_slope, linetype = 2) +
  theme_minimal() + labs(title = "Fort Severn cumulative organic carbon with depth",
    subtitle = "Two confirmed contacts; dashed pooled model is exploratory", x = "Depth (cm)", y = "kg C/m²")
ggsave(file.path(out, "cumulative_carbon_curves.png"), p, width = 8, height = 5, dpi = 180)
write_csv(data.frame(model = c("Peat depth mean", "Cumulative carbon through zero"),
  n_cores = nrow(local), estimate = c(coef(depth_model)[1], curve_slope),
  units = c("cm", "kg C/m2 per cm"),
  LOOCV_RMSE = c(sqrt(mean((depth_cv-local$peat_depth_cm)^2)),
    sqrt(mean((curve_cv-local$organic_carbon_kg_m2)^2))),
  status = "Exploratory: only two confirmed peat contacts"), file.path(out, "model_summary.csv"))
writeLines(c("Both local models ran before the Li update. No Route L models were used.",
  "Only PM-01-A and PM-03-A have explicit field-confirmed peat contacts in the inputs.",
  "All eight cores are retained in the workshop field/Sothe analysis.",
  "Add confirmed organic/mineral contact depths (including confirmed zero peat) for the other cores.",
  "Two contacts support exploratory models only; no reliable covariate depth model or validated deep carbon curve.",
  "No extrapolation of the local carbon-depth curve beyond its deepest measured organic profile.",
  "Carbon accumulation here means cumulative stock versus depth. Temporal accumulation needs chronology/ages."),
  file.path(out, "data_needs.txt"))
