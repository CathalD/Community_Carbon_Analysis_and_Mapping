# Change the input path here when substituting your own workbook.
course_config <- list(
  workbook_path = "data/input/community_cores.xlsm",
  data_origin = "Revised illustrative data",
  # Author-confirmed correction applied on import, leaving the source intact.
  density_override_J7 = 0.09,
  illustrative = TRUE,
  aoi_path = "data/aoi/study_area.geojson",
  grid_cellsize_m = 500,
  depth_datum = "surface",
  seed = 20261026
)
