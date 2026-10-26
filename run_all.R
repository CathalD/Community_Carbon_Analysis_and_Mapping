# Open the course R project, then source this file. No files outside the course are needed.
source("00_setup.R")
for (step in c("01_field", "02_compare", "03_area", "04_spatial", "05_bayes")) {
  message("Running ", step)
  source(file.path("scripts", paste0(step, ".R")))
}
capture.output(sessionInfo(), file = "outputs/R_session.txt")
message("Course complete. View docs/workshop.md and outputs/.")
