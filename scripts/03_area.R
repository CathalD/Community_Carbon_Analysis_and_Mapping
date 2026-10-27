source("R/pipeline.R")
run_stage("area", {
# Hypothetical ecosystem areas and stock totals at a shared 0–30 cm depth.
source("config.R")
source("R/spatial.R")
area_results <- run_area_analysis(course_config, cellsize_m = course_config$grid_cellsize_m)
print(area_results$area_totals, row.names = FALSE)

})
