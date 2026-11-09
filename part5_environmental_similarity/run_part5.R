# Run from the workshop root: source("part5_environmental_similarity/run_part5.R")
# Or from this folder: Rscript run_part5.R
local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- c(if (!is.null(source_file)) dirname(source_file),
                  if (length(cli)) dirname(sub("^--file=", "", cli[1])),
                  getwd(), file.path(getwd(), "part5_environmental_similarity"))
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  valid <- vapply(candidates, function(p) file.exists(file.path(p, "part5_settings.R")), logical(1))
  if (!any(valid)) stop("Cannot find the Part 5 folder.")
  bundle <- candidates[which(valid)[1]]
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(bundle)
  cat("Part 5: original-study donor comparisons. Workshop synthetic cores are separate.\n")
  analysis <- new.env(parent = globalenv())
  # Unmodified historical arithmetic-stock benchmark. Supplies reconciled sites.
  sys.source("scripts/79_compare_external_similarity_transfer.R", envir = analysis)
  sys.source("part5_settings.R", envir = analysis)
  sys.source("scripts/part5_log_donor_comparison.R", envir = analysis)
  writeLines(capture.output(sessionInfo()), "outputs/part5_log_donors/session_info.txt")
  cat("\nDone. Log donor results: ", file.path(bundle, "outputs/part5_log_donors"), "\n", sep = "")
  cat("Historical results: outputs/covariate_analysis/part15_external_similarity_transfer/\n")
  cat("Donor spread and effective donor count describe support; neither is prediction variance.\n")
})
