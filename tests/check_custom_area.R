# Portable mapping checks. All temporary inputs and outputs stay outside the course.
source("R/spatial.R")
.spatial_packages()
test_directory <- tempfile("community-carbon-mapping-")
dir.create(test_directory)

expect_mapping_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}

rectangle <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
    c(xmin, ymax), c(xmin, ymin))))
}
local_example <- function(longitude, latitude, width, height = width, multipart = FALSE) {
  crs <- paste0("+proj=laea +lat_0=", latitude, " +lon_0=", longitude,
    " +datum=WGS84 +units=m +no_defs")
  if (multipart) {
    area <- sf::st_sfc(rectangle(0, 0, width, height),
      rectangle(10 * width, 0, 11 * width, height), crs = crs)
    positions <- data.frame(x = c(.25, .75, 10.5) * width, y = c(.25, .75, .5) * height)
  } else {
    area <- sf::st_sfc(rectangle(0, 0, width, height), crs = crs)
    positions <- data.frame(x = c(.25, .5, .75) * width, y = c(.25, .75, .5) * height)
  }
  positions <- sf::st_as_sf(positions, coords = c("x", "y"), crs = crs)
  list(aoi = sf::st_sf(geometry = sf::st_transform(area, 4326)),
    coordinates = sf::st_coordinates(sf::st_transform(positions, 4326)))
}

run_mapping_case <- function(name, example) {
  folder <- file.path(test_directory, name); dir.create(folder)
  aoi_path <- file.path(folder, "my_aoi.geojson")
  sf::st_write(example$aoi, aoi_path, quiet = TRUE)
  coordinates <- example$coordinates
  cores <- data.frame(core_id = c("A", "B", "C"), ecosystem = c("Bog", "Fen", "Forest"),
    group = c("Wetland cores", "Wetland cores", "Forest cores"),
    longitude = coordinates[, 1], latitude = coordinates[, 2], stock_0_30_kg_m2 = c(6, 12, 18))
  core_path <- file.path(folder, "core_stocks.csv")
  utils::write.csv(cores, core_path, row.names = FALSE)
  config <- list(aoi_path = aoi_path, illustrative = FALSE, data_origin = "Test inputs")
  result <- run_simple_mapping(config, core_path, file.path(folder, "maps"))
  stopifnot(result$metadata$n_cores == 3, result$metadata$candidate_grid_cells <= 10000,
    nrow(result$cells) > 0, all(is.finite(result$cells$idw_stock_0_30_kg_m2)),
    all(result$cells$idw_stock_0_30_kg_m2 >= 6), all(result$cells$idw_stock_0_30_kg_m2 <= 18),
    all(file.exists(unlist(result$paths))),
    abs(sum(result$cells$area_m2) / as.numeric(sf::st_area(result$aoi)) - 1) < 1e-7)
  raster <- terra::rast(result$paths$idw_raster)
  raster_values <- terra::values(raster, mat = FALSE)
  stopifnot(terra::ncell(raster) <= 10000, any(is.finite(raster_values)),
    all(raster_values[is.finite(raster_values)] >= 6),
    all(raster_values[is.finite(raster_values)] <= 18))
  list(config = config, core_path = core_path, cores = cores, result = result, folder = folder)
}

cases <- list(tiny = local_example(-84, 53, 100),
  western = local_example(-123, 49, 50000), eastern = local_example(-63, 46, 50000),
  northern = local_example(-95, 75, 50000),
  thin = local_example(-100, 55, 50000, 1),
  multipart = local_example(-90, 55, 1000, multipart = TRUE))
for (case_name in names(cases)) checked <- run_mapping_case(case_name, cases[[case_name]])

canada_extent <- sf::st_sf(geometry = sf::st_sfc(rectangle(-140, 42, -54, 80), crs = 4326))
checked <- run_mapping_case("large_extent", list(aoi = canada_extent,
  coordinates = matrix(c(-125, 50, -90, 60, -65, 65), ncol = 2, byrow = TRUE)))

# Allocation is checked before a fine grid can consume large amounts of memory.
expect_mapping_error(.area_grid(checked$result$aoi, 1), "more than 10,000 cells")

# Invalid and incomplete cores are excluded only from mapping, with saved reasons.
checked <- run_mapping_case("exclusions", cases$western)
bad <- checked$cores
bad$core_id <- paste0(bad$core_id, "-excluded")
bad$latitude[1] <- NA_real_; bad$stock_0_30_kg_m2[2] <- NA_real_
bad$longitude[3] <- -65; bad$latitude[3] <- 45
utils::write.csv(rbind(checked$cores, bad), checked$core_path, row.names = FALSE)
result <- run_simple_mapping(checked$config, checked$core_path, file.path(checked$folder, "filtered"))
stopifnot(nrow(result$excluded) == 3, result$metadata$n_cores == 3,
  any(grepl("latitude/longitude", result$excluded$reason, fixed = TRUE)),
  any(grepl("Incomplete", result$excluded$reason)), any(grepl("Outside", result$excluded$reason)))
utils::write.csv(bad, checked$core_path, row.names = FALSE)
# Reusing a successful output folder must clear its old maps after an invalid rerun.
# A learner's unrelated file in the same folder must survive.
unrelated_file <- file.path(checked$folder, "maps", "my_notes.txt")
writeLines("Keep my notes.", unrelated_file)
expect_mapping_error(run_simple_mapping(checked$config, checked$core_path,
  file.path(checked$folder, "maps")), "at least two cores")
old_map_paths <- checked$result$paths[setdiff(names(checked$result$paths), "excluded_cores")]
current_exclusions <- utils::read.csv(checked$result$paths$excluded_cores, stringsAsFactors = FALSE)
stopifnot(!any(file.exists(unlist(old_map_paths))), nrow(current_exclusions) == 3,
  identical(current_exclusions$core_id, bad$core_id),
  identical(readLines(unrelated_file), "Keep my notes."))

duplicate <- checked$cores
duplicate$latitude[2] <- duplicate$latitude[1]; duplicate$longitude[2] <- duplicate$longitude[1]
utils::write.csv(duplicate, checked$core_path, row.names = FALSE)
expect_mapping_error(run_simple_mapping(checked$config, checked$core_path,
  file.path(checked$folder, "duplicates")), "share the same coordinates")

line_path <- file.path(checked$folder, "line.geojson")
sf::st_write(sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(rbind(c(-123, 49), c(-122, 50))),
  crs = 4326)), line_path, quiet = TRUE)
checked$config$aoi_path <- line_path
expect_mapping_error(run_simple_mapping(checked$config, checked$core_path,
  file.path(checked$folder, "line")), "polygon boundary")
message("Custom AOI checks passed: tiny, thin, multipart, western/eastern/northern Canada, large extent, allocation cap, exclusions and clear input errors.")
unlink(test_directory, recursive = TRUE)
