# Open an image already saved by a course step. No additional packages needed.
view_figure <- function(path) {
  resolved_path <- normalizePath(path, mustWork = TRUE)
  utils::browseURL(resolved_path)
  invisible(resolved_path)
}
