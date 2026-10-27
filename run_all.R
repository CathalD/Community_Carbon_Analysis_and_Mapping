# Open the course R project, then source this file. No files outside the course are needed.
source("R/pipeline.R")
invalidate_stage("field")
source("00_setup.R")
for (step in c("01_field", "02_compare", "03_area", "04_spatial", "04b_kriging", "05_bayes")) {
  message("Running ", step)
  source(file.path("scripts", paste0(step, ".R")))
}
capture.output(sessionInfo(), file = "outputs/R_session.txt")
writeLines("All course stages completed successfully.", "outputs/run_complete.txt")
message("Course complete. View docs/workshop.md and outputs/.")
