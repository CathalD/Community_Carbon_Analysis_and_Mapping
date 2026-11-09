# Optional historical multi-depth, raw-stock comparison; not a log GP.
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
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("scripts/85_run_hierarchical_similarity_fs_point_test.R", "--overwrite=true"))
  if (status != 0) stop("Hierarchical example failed; inspect the printed error.")
})
