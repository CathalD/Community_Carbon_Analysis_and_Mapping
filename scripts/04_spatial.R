source("R/pipeline.R")
run_stage("spatial", {
# IDW from four illustrative cores; separate simulated regression and validation.
source("config.R")
source("R/spatial.R")
spatial_results <- run_spatial_analysis(course_config, seed = course_config$seed)
print(spatial_results$metrics, row.names = FALSE)
print(spatial_results$totals, row.names = FALSE)

})
