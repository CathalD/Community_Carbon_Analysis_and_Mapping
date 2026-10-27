source("R/pipeline.R")
run_stage("comparison", {
# Recalculate community-to-reference comparisons and figures from the Part 1 outputs.
source("config.R")
source("R/compare.R")
comparison_results <- run_reference_comparison(course_config)

})
