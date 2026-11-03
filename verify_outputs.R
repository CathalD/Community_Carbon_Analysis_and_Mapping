library(terra)
library(readr)
source("helpers/workshop_gp.R")
checks <- list()
check <- function(label, value) {
  checks[[label]] <<- data.frame(check = label, passed = isTRUE(value))
  if (!isTRUE(value)) stop("Verification failed: ", label)
}
cores <- read_csv("outputs/script1/mapping_data.csv", show_col_types = FALSE)
layers <- read.csv("data/field/core_layers.csv")
check("All eight FS cores included", nrow(cores) == 8 && !anyDuplicated(cores[["Core ID"]]))
# Independent unit/overlap calculation from raw g/kg, g/cm3 and cm.
for (target in c(30, 50, 100)) {
  expected <- vapply(cores[["Core ID"]], function(id) {
    d <- layers[layers$plot_id == id, ]
    sum(d$soc*d$bulk_density*(1-d$coarse_frag)*pmax(0,pmin(d$depth_to,target)-d$depth_from)/100)
  }, numeric(1))
  check(paste("Measured stock integration to", target, "cm / sampled bottom"),
    max(abs(expected-cores[[paste0("C Stock 0-",target,"cm")]])) < 1e-9)
}
check("Shallow cores remain included in S30", all(is.finite(cores[["C Stock 0-30cm"]])))
raw <- readRDS("outputs/script5/sothe/stock_0_30/raw_cores/raw_core_results.rds")
check("Sothe update uses all eight cores", nrow(raw$cores) == 8)
check("Finite Sothe posterior", all(is.finite(as.matrix(raw$predictions))))
check("Nonnegative contracted Sothe log variance", all(raw$posterior_log_var >= 0) &&
  all(raw$predictions$log_variance_reduction >= -1e-8 & raw$predictions$log_variance_reduction <= 1+1e-8))
check("All 119 available covariate bands retained", nlyr(rast("data/covariates/all_predictors_1km.tif")) == 119)
# Numerical test of conditioning: one exact observation must be reproduced.
p <- workshop_gp(matrix(c(0,0),1), 5, 10, 2, matrix(c(0,0),1), 10, 2)
check("GP exact-observation conditioning", abs(p$posterior_mean-5) < 1e-10 && p$posterior_sd < 1e-8)
source("validate_sothe.R")
check("Factored Li engine reproduces workshop Sothe baseline", TRUE)
source("verify_full_suite.R")
check("Four-method Sothe and Li suites verified", all(full_suite_checks$passed))
dir.create("outputs/verification", recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, checks), "outputs/verification/checks.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "outputs/verification/session_info.txt")
code <- list.files(".", pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
code <- code[!grepl("archive/", code)]
for (f in code) parse(f)
write.csv(data.frame(file = code, md5 = unname(tools::md5sum(code))), "outputs/verification/code_manifest.csv", row.names = FALSE)
inputs <- list.files("data", full.names = TRUE, recursive = TRUE)
inputs <- inputs[!dir.exists(inputs)]
write.csv(data.frame(file = inputs, bytes = file.info(inputs)$size,
  md5 = unname(tools::md5sum(inputs))), "outputs/verification/input_manifest.csv", row.names = FALSE)
message("All FS workshop verification checks passed.")
