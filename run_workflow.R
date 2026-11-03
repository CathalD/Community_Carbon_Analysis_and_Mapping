# Open FS_test.Rproj, then source("run_workflow.R"). Reruns replace FS test outputs only.
run_fs_workflow <- function(steps = c("prepare", "1", "2", "3", "4", "5", "5c", "3d", "5d", "verify", "results")) {
  scripts <- c(prepare = "prepare_inputs.R", `1` = "script_1.R", `2` = "script_2.R",
    `3` = "script_3.R", `4` = "script_4.R", `5` = "script_5.R",
    `5c` = "script_5c_sothe_regression_update.R",
    `3d` = "script_3d_li_model_suite.R", `5d` = "script_5d_li_update_suite.R",
    `3b` = "script_3b_peat_depth_and_curves.R", `5b` = "script_5b_li_update.R",
    verify = "verify_outputs.R", results = "show_full_suite_results.R")
  stopifnot(all(steps %in% names(scripts)), file.exists("FS_test.Rproj"))
  dir.create("outputs/logs", recursive = TRUE, showWarnings = FALSE)
  status <- list()
  if (file.exists("outputs/run_status.csv")) {
    previous <- read.csv("outputs/run_status.csv", stringsAsFactors = FALSE)
    previous <- previous[previous$step %in% names(scripts) & previous$exit_code %in% c(0, 1), ]
    for (s in unique(previous$step)) status[[s]] <- tail(previous[previous$step == s, ], 1)
  }
  for (step in steps) {
    script <- scripts[[step]]
    message("FS workshop step ", step, ": ", script)
    started <- Sys.time()
    # Fresh process for every script; working directory remains FS test.
    logfile <- file.path("outputs/logs", paste0("step_", step, ".log"))
    code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
      stdout = logfile, stderr = logfile)
    status[[step]] <- data.frame(step = step, script = script, exit_code = code,
      started = as.character(started), ended = as.character(Sys.time()))
    write.csv(do.call(rbind, status), "outputs/run_status.csv", row.names = FALSE)
    if (code != 0) stop("Step ", step, " failed; inspect outputs/logs/step_", step, ".log")
  }
  invisible(do.call(rbind, status))
}
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) run_fs_workflow(args) else run_fs_workflow()
}
