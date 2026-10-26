# Simple normal–normal update with explicitly assumed, illustrative uncertainty.
source("config.R")
source("R/spatial.R")
bayesian_results <- run_bayesian_update(course_config)
print(bayesian_results$summary, row.names = FALSE)
