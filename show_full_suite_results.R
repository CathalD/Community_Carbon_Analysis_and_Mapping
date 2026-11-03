# Static figures from the completed four-method suites. No analysis is refitted.
library(terra)
library(sf)
library(ggplot2)
library(scales)
out <- "outputs/results/full_suite"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
boundary <- st_transform(st_read("data/boundary/aoi.geojson", quiet = TRUE), 3978)
field <- read.csv("outputs/script3d_li_suite/field_observations.csv")
field <- field[field$canfit, , drop = FALSE]
field$contact_evidence <- ifelse(field$confirmed, "Field-confirmed", "Layer-classified")
points <- st_transform(st_as_sf(field, coords = c("longitude", "latitude"), crs = 4326), 3978)
methods <- c("raw_cores", "area_based", "idw", "regression")
method_labels <- c(raw_cores = "Spatial core update", area_based = "Area-based update",
  idw = "IDW update", regression = "Covariate-regression update")
summaries <- list()

map_values <- function(r) {
  d <- as.data.frame(r[[1]], xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", "value")
  d
}

panel_data <- function(rasters, labels, units, count_negatives = FALSE) {
  summary <- do.call(rbind, lapply(seq_along(rasters), function(i) {
    v <- values(rasters[[i]], mat = FALSE); v <- v[is.finite(v)]
    data.frame(method = names(rasters)[i], method_label = labels[i], valid_cells = length(v),
      mean = mean(v), minimum = min(v), maximum = max(v), negative_cells = sum(v < 0))
  }))
  summary$panel <- paste0(summary$method_label, "\nMean ", formatC(summary$mean, digits = 2, format = "f"),
    " ", units, if (count_negatives) paste0("\nNegative cells: ", summary$negative_cells) else "")
  data <- do.call(rbind, lapply(seq_along(rasters), function(i) {
    d <- map_values(rasters[[i]])
    d$panel <- summary$panel[i]
    d
  }))
  data$panel <- factor(data$panel, levels = summary$panel)
  list(data = data, summary = summary, panels = summary$panel)
}

fill_scale <- function(values, units) {
  limits <- range(values, finite = TRUE)
  if (diff(limits) == 0) limits <- limits + c(-.5, .5)
  if (limits[1] < 0 && limits[2] > 0) {
    breaks <- sort(unique(c(limits[1], 0, pretty(c(0, limits[2]), n = 4))))
    breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
    scale_fill_gradientn(colours = c("#2166ac", "#f7f7f7", "#fdae61", "#a50026"),
      values = scales::rescale(c(limits[1], 0, limits[2]*.5, limits[2])), limits = limits,
      breaks = breaks, labels = scales::label_number(accuracy = .1), name = units, na.value = "white")
  } else scale_fill_viridis_c(option = "C", limits = limits, name = units, na.value = "white")
}

base_map <- function(title, subtitle, caption, ncol = 3) {
  ggplot() +
    geom_sf(data = boundary, fill = NA, colour = "#71717a", linewidth = .25) +
    coord_sf(crs = st_crs(3978), datum = NA, expand = FALSE) +
    facet_wrap(~panel, ncol = ncol) +
    scale_shape_manual(values = c("Field-confirmed" = 21, "Layer-classified" = 24),
      name = "Organic-contact evidence") +
    labs(title = title, subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      strip.text = element_text(size = 11.5, face = "bold", margin = margin(5, 1, 5, 1)),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 11, hjust = 0),
      legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 11),
      legend.text = element_text(size = 10), plot.margin = margin(12, 12, 12, 12))
}

save_plot <- function(p, filename, width = 14, height = 11) {
  ggsave(file.path(out, paste0(filename, ".png")), p, width = width, height = height,
    dpi = 160, bg = "white")
  ggsave(file.path(out, paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
}

make_update_figure <- function(rasters, labels, units, filename, title, subtitle, caption,
                               include_negatives = FALSE) {
  bundle <- panel_data(rasters, labels, units, include_negatives)
  pt <- do.call(rbind, lapply(bundle$panels, function(panel) {
    x <- points; x$panel <- factor(panel, levels = bundle$panels); x
  }))
  p <- base_map(title, subtitle, caption) +
    geom_raster(data = bundle$data, aes(x, y, fill = value)) +
    geom_sf(data = boundary, fill = NA, colour = "#71717a", linewidth = .25) +
    geom_sf(data = pt, aes(shape = contact_evidence), fill = "black", colour = "white",
      stroke = .4, size = 1.8) +
    fill_scale(bundle$data$value, units) +
    guides(fill = guide_colourbar(barwidth = grid::unit(9, "cm"), barheight = grid::unit(.4, "cm")))
  save_plot(p, filename)
  write.csv(bundle$summary, file.path(out, paste0(filename, "_summary.csv")), row.names = FALSE)
  summaries[[filename]] <<- bundle$summary
}

# Sothe prior plus four updates, all on one colour scale.
sothe_path <- "outputs/script5/sothe/stock_0_30"
sothe_spatial <- rast(file.path(sothe_path, "raw_cores/raw_core_maps.tif"))
sothe <- list(prior = sothe_spatial$prior_mean, raw_cores = sothe_spatial$posterior_mean,
  area_based = rast(file.path(sothe_path, "area_based/posterior_maps.tif"))$posterior_mean,
  idw = rast(file.path(sothe_path, "idw/posterior_maps.tif"))$posterior_mean,
  regression = rast(file.path(sothe_path, "regression/posterior_maps.tif"))$posterior_mean)
make_update_figure(sothe, c("Sothe prior", unname(method_labels)), "kg C/m²", "sothe_prior_and_four_updates",
  "Fort Severn | Sothe prior and all four updates",
  paste0(nrow(field), " FS cores included; one shared stock scale; 1 km analysis grid."),
  "Circles: field-confirmed organic contacts; triangles: layer-classified contacts. All eight enter Sothe stock calculations.\nThese are alternative estimates, not measured carbon loss. S30 uncertainty uses the requested S100-derived proxy.")

for (variable in c("depth", "carbon")) {
  units <- if (variable == "depth") "cm" else "kg C/m²"
  quantity <- if (variable == "depth") "peat depth" else "organic carbon"
  li_path <- file.path("outputs/script5/li_suite", variable)
  updates <- lapply(methods, function(method) rast(file.path(li_path, method, "posterior_maps.tif")))
  names(updates) <- methods
  li <- c(list(prior = updates$raw_cores$prior_mean), lapply(updates, function(r) r$posterior_mean))
  make_update_figure(li, c("Li prior", unname(method_labels)), units,
    paste0("li_", variable, "_prior_and_four_updates"),
    paste0("Fort Severn | Li ", quantity, ": prior and four updates"),
    "Original-scale Gaussian sensitivity: all eight contact labels, including four exact zeros.",
    "Circles: confirmed contacts; triangles: layer-classified contacts. Negative estimates are blue, retained and counted.\nNegative values are model diagnostics, not physical peat depth/carbon. Regression retains the prior where local estimates are negative/missing.\nDepth and carbon updates are dependent comparisons; they are not a joint inventory uncertainty.",
    include_negatives = TRUE)

  # Four local methods: observations are shown only as points, never a fake raster.
  local <- rast(file.path("outputs/script3d_li_suite", variable, "local_maps.tif"))
  local_rasters <- lapply(c("area_based", "idw", "regression"), function(method) local[[method]])
  names(local_rasters) <- c("area_based", "idw", "regression")
  labels <- c("Area-average model", "IDW model", "Covariate-regression model")
  bundle <- panel_data(local_rasters, labels, units, count_negatives = TRUE)
  raw_label <- paste0("Raw core observations\nPoint mean ", formatC(mean(field[[variable]]), digits = 2, format = "f"),
    " ", units, "\nn = ", nrow(field))
  panels <- c(raw_label, bundle$panels)
  bundle$data$panel <- factor(as.character(bundle$data$panel), levels = panels)
  raw_points <- points
  raw_points$observed <- field[[variable]]
  raw_points$panel <- factor(raw_label, levels = panels)
  other_points <- do.call(rbind, lapply(bundle$panels, function(panel) {
    x <- points; x$panel <- factor(panel, levels = panels); x
  }))
  p <- base_map(paste0("Fort Severn | Li ", quantity, ": four local methods"),
    "Raw measured/classified points and three local maps, before using Li priors.",
    "The raw panel contains points only. Circles: confirmed contacts; triangles: layer-classified contacts.\nOne shared scale includes negative extrapolations (blue), retained and flagged.\nMean labels describe points in the raw panel and mapped cells in the other panels.",
    ncol = 2) +
    geom_raster(data = bundle$data, aes(x, y, fill = value)) +
    geom_sf(data = boundary, fill = NA, colour = "#71717a", linewidth = .25) +
    geom_sf(data = other_points, aes(shape = contact_evidence), fill = "black", colour = "white", stroke = .4, size = 1.6) +
    geom_sf(data = raw_points, aes(shape = contact_evidence, fill = observed), colour = "black", stroke = .5, size = 3) +
    fill_scale(c(bundle$data$value, field[[variable]]), units) +
    guides(fill = guide_colourbar(barwidth = grid::unit(9, "cm"), barheight = grid::unit(.4, "cm")))
  save_plot(p, paste0("li_", variable, "_four_local_methods"), width = 11, height = 11)
  raw_summary <- data.frame(method = "raw_cores", method_label = "Raw core observations",
    valid_cells = NA_integer_, mean = mean(field[[variable]]), minimum = min(field[[variable]]),
    maximum = max(field[[variable]]), negative_cells = sum(field[[variable]] < 0), panel = raw_label)
  write.csv(rbind(raw_summary, bundle$summary),
    file.path(out, paste0("li_", variable, "_four_local_methods_summary.csv")), row.names = FALSE)
}

file.copy("outputs/script3d_li_suite/carbon_curve/cumulative_carbon_curves.png",
  file.path(out, "li_cumulative_carbon_curves.png"), overwrite = TRUE)
write.csv(data.frame(figure = list.files(out, pattern = "[.]png$"),
  description = c("Li carbon: four local methods", "Li carbon: prior plus four updates",
    "Local cumulative carbon curves", "Li depth: four local methods",
    "Li depth: prior plus four updates", "Sothe: prior plus four updates")),
  file.path(out, "figure_index.csv"), row.names = FALSE)
message("Full four-method presentation figures generated in ", out)
