required_packages <- c("readxl", "ggplot2", "sf", "terra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
if (!all(vapply(required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1)))) stop("Package setup is incomplete. Review installation messages.")
if (!file.exists("config.R")) stop("Open the course R project before running this script.")
message("Setup complete: ", paste(required_packages, collapse = ", "))
