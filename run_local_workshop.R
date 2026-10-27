# Start with field + AOI. Enable each later stage only after preparing its inputs.
source("00_setup.R")
source("R/local.R")
local_settings <- list(
  field_path = "data/input/my_field_data.csv",
  aoi_path = "data/aoi/my_aoi.geojson",
  output_dir = "outputs/my_data",
  cellsize_m = NULL,
  ecosystem_cells = NULL,       # e.g. data/input/my_ecosystem_cells.csv
  predictor_cells = NULL,       # e.g. data/input/my_predictor_cells.csv
  model_observations = NULL,    # e.g. data/input/my_model_observations.csv
  holdout_buffer_m = NULL,      # choose from spatial dependence and intended prediction task
  prior_cells = NULL,           # e.g. data/input/my_prior_cells.csv
  independent_errors_confirmed = FALSE
)
local_results <- run_local_workshop(local_settings)
