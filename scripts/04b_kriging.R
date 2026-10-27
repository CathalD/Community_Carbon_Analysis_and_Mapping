source("R/pipeline.R")
run_stage("kriging", {
source("config.R")
source("R/kriging.R")
kriging_metrics <- run_kriging_demo(course_config)
print(kriging_metrics)

})
