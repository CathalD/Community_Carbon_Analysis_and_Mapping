# Start here. The included example needs no file changes.
source("00_setup.R")
source("R/course.R")
example_results <- run_course(
  workbook_path = "data/input/community_cores.xlsm",
  aoi_path = "data/aoi/study_area.geojson",
  output_dir = "outputs/example",
  density_override_J7 = 0.09,
  data_origin = "Revised illustrative data"
)
