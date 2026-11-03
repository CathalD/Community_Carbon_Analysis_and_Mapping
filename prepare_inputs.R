# Make small, local FS inputs from existing rasters without copying regional files.
library(terra)
source("config.R")
terra::projNetwork(FALSE)
dir.create("outputs/preflight", recursive = TRUE, showWarnings = FALSE)
aoi <- project(vect("data/boundary/aoi.geojson"), fs_config$analysis_crs)
step <- fs_config$analysis_resolution_m
template <- rast(xmin = floor(xmin(aoi)/step)*step, xmax = ceiling(xmax(aoi)/step)*step,
  ymin = floor(ymin(aoi)/step)*step, ymax = ceiling(ymax(aoi)/step)*step,
  resolution = step, crs = crs(aoi))
values(template) <- 1
template <- mask(template, aoi, touches = FALSE)
writeRaster(template, "data/boundary/analysis_grid.tif", overwrite = TRUE)
manifest <- list()
prepare_prior <- function(source, output, role) {
  if (!file.exists(source)) stop("Missing prior input: ", source)
  x <- rast(source)
  # Crop in the native grid first. Preserve each mean/spread pair with nearest transfer.
  x <- crop(x, project(aoi, crs(x)), snap = "out")
  writeRaster(x, output, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
  manifest[[length(manifest) + 1L]] <<- data.frame(role = role,
    source = normalizePath(source), local = output, bytes = file.info(source)$size,
    native_resolution = paste(res(x), collapse = " x "),
    local_md5 = unname(tools::md5sum(output)))
}
sothe <- c(sothe_s30_mean = "McMaster_WWFCanada_soil_carbon30cm_250m_kg-m2_version3.0.tif",
  sothe_s100_mean = "McMaster_WWFCanada_soil_carbon1m_250m_kg-m2_version3.0.tif",
  sothe_s100_interval90_width = "McMaster_WWFCanada_soil_carbon1m_uncertainty_250m_kg-m2_version3.0.tif")
for (nm in names(sothe)) prepare_prior(file.path(fs_config$sothe_root, sothe[[nm]]),
  file.path("data/priors", paste0(nm, ".tif")), nm)
li <- c(li_depth_mean = "McMaster_PeatDepth_30m_cm_version1.0.tif",
  li_depth_sd = "McMaster_Uncertainty_PeatDepth_30m_cm_version1.0.tif",
  li_carbon_mean = "McMaster_Carbon_30m_kgm2_version1.0.tif",
  li_carbon_sd = "McMaster_Uncertainty_Carbon_30m_kgm2_version1.0.tif")
for (nm in names(li)) prepare_prior(file.path(fs_config$source_root, "data/external", li[[nm]]),
  file.path("data/priors", paste0(nm, ".tif")), nm)
write.csv(do.call(rbind, manifest), "outputs/preflight/prior_manifest.csv", row.names = FALSE)

# All completed downloaded tiles; never use inspection samples or incomplete downloads.
# The 1 km grid matches Script 5. Covariates are bilinearly sampled at grid centres,
# not claimed to be area averages. Native values are extracted separately at cores.
locations <- read.csv("data/field/core_locations.csv")
pts <- vect(locations, geom = c("longitude", "latitude"), crs = "EPSG:4326")
predictors <- list(); core_predictors <- list(); bands <- list()
for (group in c("terrain", "surface_context", "embeddings")) {
  message("Preparing FS covariates: ", group)
  files <- list.files(fs_config$covariate_root,
    pattern = paste0("_AOI_", group, "[0-9-]+[.]tif$"), full.names = TRUE)
  if (!length(files)) stop("Missing completed covariate tiles: ", group)
  tiles <- lapply(files, rast)
  stopifnot(all(vapply(tiles, function(x) identical(names(x), names(tiles[[1]])), logical(1))))
  x <- vrt(normalizePath(files), filename = file.path("data/covariates", paste0(group, ".vrt")), overwrite = TRUE)
  names(x) <- names(tiles[[1]])
  native_values <- extract(x, project(pts, crs(x)), ID = FALSE)
  core_predictors[[group]] <- native_values
  aligned <- mask(project(x, template, method = "bilinear"), template)
  predictors[[group]] <- writeRaster(aligned, file.path("data/covariates", paste0(group, "_1km.tif")),
    overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
  bands[[group]] <- data.frame(group = group, band = names(x),
    source_tiles = paste(normalizePath(files), collapse = " | "),
    mapping_transfer = "Bilinear at 1 km grid centres; native sampling at field cores")
}
all_covariates <- do.call(c, unname(predictors))
stopifnot(!anyDuplicated(names(all_covariates)))
writeRaster(all_covariates, "data/covariates/all_predictors_1km.tif", overwrite = TRUE,
  gdal = "COMPRESS=DEFLATE")
write.csv(data.frame(core_id = locations$plot_id, do.call(cbind, unname(core_predictors))),
  "data/covariates/native_core_covariates.csv", row.names = FALSE)
write.csv(do.call(rbind, bands), "data/covariates/band_manifest.csv", row.names = FALSE)
coverage <- global(all_covariates, "notNA")
write.csv(data.frame(band = rownames(coverage), valid_cells = coverage$notNA,
  aoi_grid_cells = global(template, "notNA")[1,1]), "outputs/preflight/covariate_coverage.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "outputs/preflight/session_info.txt")
