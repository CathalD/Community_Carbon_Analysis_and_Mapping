# Put your two files at these paths. Edit only these lines for other file names.
my_field_data <- "data/input/my_field_data.csv"
my_area <- "data/aoi/my_aoi.geojson"

source("00_setup.R")
source("R/course.R")
my_results <- run_course(
  workbook_path = my_field_data,
  aoi_path = my_area,
  output_dir = "outputs/my_data"
)
