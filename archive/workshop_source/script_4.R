#===============================================================

# -------- Part 4 ------------------

#===============================================================

library(terra)
library(readr)
library(dplyr)

output_folder <- "outputs/script4"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

#============
# Step 1 - Load Everything we've done so dar
#===================

# Field stocks and study boundary

cores <- read_csv(
  "outputs/script1/mapping_data.csv",
  show_col_types = FALSE
)

aoi <- vect("data/boundary/aoi.geojson")

core_points <- vect(
  cores[, c("Core ID", "Longitude", "Latitude")],
  geom = c("Longitude", "Latitude"),
  crs = "EPSG:4326"
)


# Area-based average and IDW outputs

area_0_30 <- vect(
  "outputs/script2/area_based.gpkg",
  layer = "area_based_0_30"
)

area_0_100 <- vect(
  "outputs/script2/area_based.gpkg",
  layer = "area_based_0_100"
)

idw_0_30 <- rast("outputs/script2/idw_0_30cm.tif")
idw_0_100 <- rast("outputs/script2/idw_0_100cm.tif")

# Match the boundary and points to the IDW coordinate system
aoi <- project(aoi, crs(idw_0_30))
core_points <- project(core_points, crs(idw_0_30))

#---------------------------------------------------------------
# Script 3 regression maps — currently 0–30 cm
#---------------------------------------------------------------
# Optional selected model from Script 3
comparison_stock <- "stock_0_30"

selected_folder <- file.path(
  "outputs/script3",
  comparison_stock
)

selected_map_file <- file.path(
  selected_folder,
  "selected_carbon_stock.tif"
)

selected_info_file <- file.path(
  selected_folder,
  "selected_model.csv"
)

regression_maps <- list()

if (file.exists(selected_map_file) &&
    file.exists(selected_info_file)) {
  
  selected_info <- read_csv(
    selected_info_file,
    show_col_types = FALSE
  )
  
  model_label <- paste0("Script3_", selected_info$model[1])
  
  regression_maps[[model_label]] <- rast(selected_map_file)
  
  cat("Imported:", model_label, "\n")
  
} else {
  message("No selected Script 3 model found; continuing without it.")
}

#---------------------------------------------------------------
# Prior stock and uncertainty — 0–100 cm
#---------------------------------------------------------------

prior_0_100 <- rast(
  "data/priors/mcmaster_0_100cm.tif"
)

prior_uncertainty_0_100 <- rast(
  "data/priors/mcmaster_uncertainty_0_100cm.tif"
)

names(prior_0_100) <- "prior_stock_0_100"
names(prior_uncertainty_0_100) <- "prior_uncertainty_0_100"

print(prior_0_100)
print(prior_uncertainty_0_100)

#---------------------------------------------------------------
# Preview the prior maps with field-core locations
#---------------------------------------------------------------

prior_aoi <- project(aoi, crs(prior_0_100))
prior_points <- project(core_points, crs(prior_0_100))

old_par <- par(mfrow = c(1, 2))

plot(prior_0_100, main = "Prior stock: 0–100 cm (kg C/m²)")
plot(prior_aoi, add = TRUE, border = "black")
points(prior_points, pch = 21, bg = "red", cex = 1)

plot(prior_uncertainty_0_100, main = "Prior uncertainty")
plot(prior_aoi, add = TRUE, border = "black")
points(prior_points, pch = 21, bg = "red", cex = 1)

par(old_par)


#==========================================
#----- 4B: Align maps for comparison -------
#==========================================

# Convert the area-average polygons to rasters
area_raster_0_30 <- rasterize(
  project(area_0_30, crs(idw_0_30)),
  idw_0_30,
  field = "mean_stock_kg_m2"
)

area_raster_0_100 <- rasterize(
  project(area_0_100, crs(idw_0_100)),
  idw_0_100,
  field = "mean_stock_kg_m2"
)


# Align a stock map to the comparison grid
align_stock <- function(x, template) {
  
  if (compareGeom(x, template, stopOnError = FALSE)) {
    return(x)
  }
  
  if (crs(x) != crs(template)) {
    return(project(x, template, method = "average"))
  }
  
  resample(x, template, method = "average")
}


# Combine maps and retain their shared coverage
prepare_comparison <- function(map_list, template, boundary) {
  
  aligned <- lapply(
    map_list,
    function(x) align_stock(x, template)
  )
  
  maps <- do.call(c, unname(aligned))
  names(maps) <- names(map_list)
  
  maps <- mask(
    maps,
    project(boundary, crs(template)),
    touches = FALSE
  )
  
  # A cell qualifies only when all methods have data
  shared_cells <- app(!is.na(maps), fun = "sum") == nlyr(maps)
  
  if (global(shared_cells, "sum", na.rm = TRUE)[1, 1] == 0) {
    stop("The maps have no shared coverage inside the boundary.")
  }
  
  mask(maps, shared_cells, maskvalues = 0)
}


#----- 0–30 cm comparison -----

maps_0_30 <- c(
  list(
    Area_average = area_raster_0_30,
    IDW = idw_0_30
  ),
  regression_maps
)

comparison_0_30 <- prepare_comparison(
  maps_0_30,
  idw_0_30,
  aoi
)


#----- 0–100 cm comparison -----

maps_0_100 <- list(
  Area_average = area_raster_0_100,
  IDW = idw_0_100,
  McMaster_prior = prior_0_100
)

comparison_0_100 <- prepare_comparison(
  maps_0_100,
  idw_0_100,
  aoi
)


#----- Save the aligned maps -----

writeRaster(
  comparison_0_30,
  file.path(output_folder, "comparison_0_30.tif"),
  overwrite = TRUE
)

writeRaster(
  comparison_0_100,
  file.path(output_folder, "comparison_0_100.tif"),
  overwrite = TRUE
)


#----- Check the number of comparison cells -----

cat("\n0–30 cm: cells available for each method\n")
print(global(!is.na(comparison_0_30), "sum"))

cat("\n0–100 cm: cells available for each method\n")
print(global(!is.na(comparison_0_100), "sum"))


#==========================================
#----- 4C: Summary statistics and chart ----
#==========================================

library(ggplot2)

comparison_sets <- list(
  "0-30 cm" = comparison_0_30,
  "0-100 cm" = comparison_0_100
)

depth_files <- c(
  "0-30 cm" = "0_30",
  "0-100 cm" = "0_100"
)

# All stock statistics are in kg C/m².
# n_cells counts map cells, not field observations.

map_summary <- data.frame()

for (depth in names(comparison_sets)) {
  
  r <- comparison_sets[[depth]]
  
  for (method in names(r)) {
    
    x <- values(r[[method]], mat = FALSE)
    x <- x[is.finite(x)]
    
    map_summary <- bind_rows(
      map_summary,
      data.frame(
        depth = depth,
        method = method,
        n_cells = length(x),
        mean_stock = mean(x),
        median_stock = median(x),
        spatial_sd = sd(x),
        spatial_cv_percent = if (mean(x) == 0) {
          NA_real_
        } else {
          100 * sd(x) / abs(mean(x))
        },
        q25 = unname(quantile(x, 0.25)),
        q75 = unname(quantile(x, 0.75)),
        minimum = min(x),
        maximum = max(x)
      )
    )
  }
}

print(map_summary)

write_csv(
  map_summary,
  file.path(output_folder, "map_summary.csv")
)


# Horizontal chart: mean and middle 50% of map values

summary_plot <- ggplot(
  map_summary,
  aes(y = method)
) +
  geom_segment(
    aes(x = q25, xend = q75, yend = method),
    linewidth = 3,
    colour = "#A6CEE3"
  ) +
  geom_point(
    aes(x = mean_stock),
    size = 3,
    colour = "#1F78B4"
  ) +
  facet_wrap(~ depth, ncol = 1, scales = "free") +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Carbon stock across the comparison area",
    x = "Carbon stock (kg C/m²)",
    y = NULL,
    caption = paste(
      "Points: map means. Bars: 25th–75th percentiles.",
      "Bars describe spatial variation, not confidence intervals."
    )
  ) +
  theme_minimal(base_size = 12)

print(summary_plot)

ggsave(
  file.path(output_folder, "map_summary.png"),
  summary_plot,
  width = 9,
  height = 6,
  dpi = 300
)


#==========================================
#----- 4D: Map comparison panels -----------
#==========================================

# Each panel uses the same colour scale within its depth group.
# Black points mark field cores.

save_map_panel <- function(
    r, depth, filename, difference = FALSE) {
  
  limits <- range(minmax(r), finite = TRUE)
  
  if (difference) {
    
    largest <- max(abs(limits))
    
    if (largest == 0) largest <- 1
    
    limits <- c(-largest, largest)
    
    colours <- colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(101)
    
  } else {
    
    if (diff(limits) == 0) {
      limits <- limits + c(-0.5, 0.5)
    }
    
    colours <- hcl.colors(100, "YlGnBu")
  }
  
  boundary <- project(aoi, crs(r))
  locations <- project(core_points, crs(r))
  
  png(
    file.path(output_folder, filename),
    width = 650 * nlyr(r),
    height = 700,
    res = 150
  )
  
  on.exit(dev.off())
  
  par(mfrow = c(1, nlyr(r)))
  
  for (i in seq_len(nlyr(r))) {
    
    label <- paste(
      strwrap(gsub("_", " ", names(r)[i]), width = 25),
      collapse = "\n"
    )
    
    plot(
      r[[i]],
      main = paste(depth, label, sep = "\n"),
      range = limits,
      col = colours,
      plg = list(title = "kg C/m²"),
      mar = c(2, 2, 5, 5)
    )
    
    plot(
      boundary,
      add = TRUE,
      border = "grey30",
      lwd = 0.7
    )
    
    points(
      locations,
      pch = 21,
      bg = "black",
      col = "white",
      cex = 1
    )
  }
}

for (depth in names(comparison_sets)) {
  
  save_map_panel(
    comparison_sets[[depth]],
    depth,
    paste0("maps_", depth_files[[depth]], ".png")
  )
}


#==========================================
#----- 4E: Differences between methods -----
#==========================================

# Differences describe map agreement, not prediction accuracy.
# Direction is always: second method minus first method.

map_agreement <- data.frame()

for (depth in names(comparison_sets)) {
  
  r <- comparison_sets[[depth]]
  pairs <- combn(names(r), 2)
  difference_layers <- list()
  
  for (i in seq_len(ncol(pairs))) {
    
    first <- pairs[1, i]
    second <- pairs[2, i]
    
    a <- values(r[[first]], mat = FALSE)
    b <- values(r[[second]], mat = FALSE)
    
    keep <- is.finite(a) & is.finite(b)
    difference_values <- b[keep] - a[keep]
    
    label <- paste0(second, "_minus_", first)
    
    difference_layers[[label]] <- r[[second]] - r[[first]]
    
    map_agreement <- bind_rows(
      map_agreement,
      data.frame(
        depth = depth,
        first_method = first,
        second_method = second,
        n_cells = sum(keep),
        mean_difference = mean(difference_values),
        mean_absolute_difference = mean(abs(difference_values)),
        rms_difference = sqrt(mean(difference_values^2))
      )
    )
  }
  
  difference_maps <- do.call(c, unname(difference_layers))
  names(difference_maps) <- names(difference_layers)
  
  writeRaster(
    difference_maps,
    file.path(
      output_folder,
      paste0("differences_", depth_files[[depth]], ".tif")
    ),
    overwrite = TRUE
  )
  
  save_map_panel(
    difference_maps,
    depth,
    paste0("differences_", depth_files[[depth]], ".png"),
    difference = TRUE
  )
}

print(map_agreement)

write_csv(
  map_agreement,
  file.path(output_folder, "map_agreement.csv")
)


#==========================================
#----- 4F: Compare maps at field cores -----
#==========================================

# Extract the mapped cell containing each core.
# Keep the field stock definitions established in Script 1.

field_columns <- c(
  "0-30 cm" = "C Stock 0-30cm",
  "0-100 cm" = "C Stock 0-100cm"
)

core_comparison <- data.frame()

for (depth in names(comparison_sets)) {
  
  r <- comparison_sets[[depth]]
  
  mapped_values <- terra::extract(
    r,
    project(core_points, crs(r)),
    method = "simple",
    ID = FALSE
  )
  
  for (method in names(r)) {
    
    core_comparison <- bind_rows(
      core_comparison,
      data.frame(
        core_id = cores[["Core ID"]],
        depth = depth,
        method = method,
        observed = cores[[field_columns[[depth]]]],
        mapped = mapped_values[[method]]
      )
    )
  }
}

core_comparison <- core_comparison %>%
  mutate(difference = mapped - observed)

write_csv(
  core_comparison,
  file.path(output_folder, "core_comparison.csv")
)


# Descriptive fit statistics at the sampled cores

core_fit_summary <- core_comparison %>%
  filter(is.finite(observed), is.finite(mapped)) %>%
  group_by(depth, method) %>%
  summarise(
    n_cores = n(),
    mean_observed = mean(observed),
    mean_mapped = mean(mapped),
    mean_difference = mean(difference),
    MAE = mean(abs(difference)),
    RMSE = sqrt(mean(difference^2)),
    .groups = "drop"
  )

print(core_fit_summary)

write_csv(
  core_fit_summary,
  file.path(output_folder, "core_fit_summary.csv")
)


# Observed versus mapped stocks

core_plot <- ggplot(
  core_comparison,
  aes(x = observed, y = mapped)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    size = 2.5,
    colour = "#1F78B4",
    na.rm = TRUE
  ) +
  facet_wrap(~ depth + method, ncol = 3) +
  coord_equal() +
  labs(
    title = "Mapped carbon stock at field-core locations",
    x = "Observed stock (kg C/m²)",
    y = "Mapped stock (kg C/m²)",
    caption = paste(
      "Dashed line: mapped equals observed.",
      "Local-method fits are not independent validation."
    )
  ) +
  theme_minimal(base_size = 11)

print(core_plot)

ggsave(
  file.path(output_folder, "core_comparison.png"),
  core_plot,
  width = 11,
  height = 7,
  dpi = 300
)


#==========================================
#----- 4G: Carry prior files into Script 5 -
#==========================================

# Copy both prior files unchanged.
# The uncertainty values are not converted or reinterpreted here.

prior_files <- c(
  "mcmaster_0_100cm.tif",
  "mcmaster_uncertainty_0_100cm.tif"
)

copied <- file.copy(
  from = file.path("data/priors", prior_files),
  to = output_folder,
  overwrite = TRUE
)

stopifnot(all(copied))

cat(
  "\nScript 4 complete. Outputs saved in:\n",
  normalizePath(output_folder),
  "\n"
)

