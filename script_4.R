# Workshop Part 4: align and compare the FS area / IDW / covariate maps and Sothe prior.
library(terra)
library(readr)
library(dplyr)
source("config.R")
dir.create("outputs/script4", recursive = TRUE, showWarnings = FALSE)
template <- rast("data/boundary/analysis_grid.tif")
aoi <- project(vect("data/boundary/aoi.geojson"), crs(template))
align <- function(x, method = "near") mask(project(x, template, method = method), template)
prior <- rast("data/priors/sothe_s30_mean.tif")
width <- rast("data/priors/sothe_s100_interval90_width.tif")
names(prior) <- "Sothe_S30_prior"
names(width) <- "Sothe_S100_full_90pct_width_proxy"
file.copy("data/priors/sothe_s30_mean.tif", "outputs/script4/sothe_s30_mean.tif", overwrite = TRUE)
file.copy("data/priors/sothe_s100_interval90_width.tif", "outputs/script4/sothe_s100_interval90_width.tif", overwrite = TRUE)
area <- project(vect("outputs/script2/area_based.gpkg", layer = "area_based_0_30"), crs(template))
area_raster <- rasterize(area, template, field = "mean_stock_kg_m2")
idw <- align(rast("outputs/script2/idw_0_30cm.tif"), "bilinear")
maps <- c(area_raster, idw, align(prior))
names(maps) <- c("Area_average", "IDW", "Sothe_S30_prior")
selected_file <- "outputs/script3/stock_0_30/selected_carbon_stock.tif"
if (file.exists(selected_file)) {
  selected <- align(rast(selected_file), "bilinear")
  names(selected) <- "Script3_selected"
  maps <- c(maps, selected)
}
common <- app(is.finite(maps), "sum") == nlyr(maps)
stopifnot(global(common, "sum", na.rm = TRUE)[1,1] > 0)
maps <- mask(maps, common, maskvalues = 0)
writeRaster(maps, "outputs/script4/comparison_0_30.tif", overwrite = TRUE)
writeRaster(align(width), "outputs/script4/s100_width_proxy_aligned.tif", overwrite = TRUE)
summary <- bind_rows(lapply(names(maps), function(nm) {
  v <- values(maps[[nm]], mat = FALSE); v <- v[is.finite(v)]
  data.frame(depth = "0-30 cm / sampled bottom", method = nm, n_cells = length(v),
    mean_stock = mean(v), spatial_sd = sd(v), minimum = min(v), maximum = max(v))
}))
write_csv(summary, "outputs/script4/map_summary.csv")
cores <- read_csv("outputs/script1/mapping_data.csv", show_col_types = FALSE)
pts <- project(vect(cores, geom = c("Longitude", "Latitude"), crs = "EPSG:4326"), crs(maps))
at_core <- extract(maps, pts, ID = FALSE)
write_csv(bind_cols(cores %>% select(`Core ID`, `C Stock 0-30cm`), at_core), "outputs/script4/core_comparison.csv")
pairs <- combn(names(maps), 2, simplify = FALSE)
differences <- lapply(pairs, function(p) maps[[p[2]]] - maps[[p[1]]])
diff <- do.call(c, differences)
names(diff) <- vapply(pairs, function(p) paste(p[2], "minus", p[1], sep = "_"), character(1))
writeRaster(diff, "outputs/script4/differences_0_30.tif", overwrite = TRUE)
png("outputs/script4/maps_0_30.png", width = 1800, height = 1400, res = 160)
plot(maps, range = range(minmax(maps), finite = TRUE)); dev.off()
writeLines(c("Sothe S30 mean is the 0-30 cm stock prior.",
  "Sothe S100 uncertainty is a full Q95-Q05 interval width, not SD.",
  "Script 5 uses width/(2*qnorm(.95)) as an assumed absolute S30 SD proxy.",
  "This requested cross-depth sensitivity is not calibrated S30 uncertainty.",
  "Partial field cores are included, integrated only to their sampled bottom.",
  "Map differences and in-sample comparisons are not independent validation."),
  "outputs/script4/uncertainty_definition.txt")
