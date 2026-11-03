# Display saved results only. Does not refit or modify the analysis products.
library(terra)
library(ggplot2)
library(patchwork)
library(sf)
out <- "outputs/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 14))
boundary <- st_transform(st_read("data/boundary/aoi.geojson", quiet = TRUE), 3978)
cores <- read.csv("outputs/script1/mapping_data.csv", check.names = FALSE)
points <- st_transform(st_as_sf(cores, coords = c("Longitude", "Latitude"), crs = 4326), 3978)
contacts <- read.csv("outputs/script3b_li_local/li_field_observations.csv")
contact_points <- st_transform(st_as_sf(contacts, coords = c("longitude", "latitude"), crs = 4326), 3978)

map_plot <- function(r, title, limits, points_layer = points,
                     units = "kg C/m²", difference = FALSE, subtitle = NULL) {
  d <- as.data.frame(r[[1]], xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", "value")
  p <- ggplot(d, aes(x, y, fill = value)) + geom_raster() +
    geom_sf(data = boundary, inherit.aes = FALSE, fill = NA, colour = "#71717a", linewidth = .25) +
    geom_sf(data = points_layer, inherit.aes = FALSE, shape = 21, fill = "black",
      colour = "white", stroke = .4, size = 1.9) +
    coord_sf(crs = st_crs(3978), datum = NA, expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, fill = units) +
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      plot.title = element_text(size = 15, face = "bold"), plot.subtitle = element_text(size = 12),
      legend.position = "bottom", legend.key.width = grid::unit(.8, "cm"),
      legend.title = element_text(size = 11), legend.text = element_text(size = 10),
      plot.margin = margin(6, 8, 6, 8))
  if (difference) p + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = limits)
  else p + scale_fill_viridis_c(option = "C", limits = limits)
}

sothe_path <- "outputs/script5/sothe/stock_0_30"
spatial <- rast(file.path(sothe_path, "raw_cores/raw_core_maps.tif"))
area <- rast(file.path(sothe_path, "area_based/posterior_maps.tif"))
idw <- rast(file.path(sothe_path, "idw/posterior_maps.tif"))
limits <- c(0, max(global(spatial$prior_mean, "max", na.rm = TRUE)[1,1],
  global(spatial$posterior_mean, "max", na.rm = TRUE)[1,1],
  global(area$posterior_mean, "max", na.rm = TRUE)[1,1],
  global(idw$posterior_mean, "max", na.rm = TRUE)[1,1]))
sothe <- (map_plot(spatial$prior_mean, "Sothe prior", limits, subtitle = "Map mean: 31.39 kg C/m²") +
  map_plot(spatial$posterior_mean, "Spatial core update", limits, subtitle = "Map mean: 29.98 kg C/m²")) /
  (map_plot(area$posterior_mean, "Area-based update", limits, subtitle = "Map mean: 13.18 kg C/m²") +
  map_plot(idw$posterior_mean, "IDW update", limits, subtitle = "Map mean: 13.99 kg C/m²")) +
  plot_annotation(title = "Fort Severn | Sothe 0–30 cm workshop results",
    subtitle = "All eight cores included to the target or sampled bottom; the same colour scale in every map.",
    caption = "Black dots: core locations. White gaps: NoData. Map means use the same 3,873 cells.\nThese are alternative estimates, not measured carbon loss. Uncertainty uses an S100-derived S30 proxy.",
    theme = theme(plot.title = element_text(size = 19, face = "bold"),
      plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 11, hjust = 0)))
ggsave(file.path(out, "sothe_map_comparison.png"), sothe, width = 11, height = 12, dpi = 160, bg = "white")
ggsave(file.path(out, "sothe_map_comparison.pdf"), sothe, width = 11, height = 12, device = cairo_pdf)

# Separate change from conditional uncertainty to avoid interpreting one as the other.
delta_limit <- max(abs(values(spatial$change_from_prior)), na.rm = TRUE)
sd_limit <- max(values(spatial$prior_sd), values(spatial$posterior_sd), na.rm = TRUE)
unc <- map_plot(spatial$change_from_prior, "Spatial update minus prior", c(-delta_limit, delta_limit), difference = TRUE) +
  map_plot(spatial$prior_sd, "Prior SD proxy", c(0, sd_limit)) +
  map_plot(spatial$posterior_sd, "Updated SD proxy", c(0, sd_limit)) +
  plot_annotation(title = "Sothe spatial update | change and uncertainty proxy",
    caption = "The uncertainty panels are conditional on transferring the S100 interval width to S30; they are not calibrated S30 uncertainty.",
    theme = theme(plot.title = element_text(size = 18, face = "bold"), plot.caption = element_text(hjust = 0)))
ggsave(file.path(out, "sothe_change_and_uncertainty.png"), unc, width = 15, height = 6, dpi = 160, bg = "white")

depth <- rast("outputs/script5/li/depth/posterior_maps.tif")
carbon <- rast("outputs/script5/li/carbon/posterior_maps.tif")
dlim <- c(0, max(values(depth$prior_mean), values(depth$posterior_mean), na.rm = TRUE))
clim <- c(0, max(values(carbon$prior_mean), values(carbon$posterior_mean), na.rm = TRUE))
li <- (map_plot(depth$prior_mean, "Li peat-depth prior", dlim, contact_points, "cm", subtitle = "Map mean: 125.33 cm") +
  map_plot(depth$posterior_mean, "Peat-depth update", dlim, contact_points, "cm", subtitle = "Map mean: 124.83 cm")) /
  (map_plot(carbon$prior_mean, "Li organic-carbon prior", clim, contact_points, subtitle = "Map mean: 64.14 kg C/m²") +
  map_plot(carbon$posterior_mean, "Organic-carbon update", clim, contact_points, subtitle = "Map mean: 63.50 kg C/m²")) +
  plot_annotation(title = "Fort Severn | Li updates — exploratory",
    subtitle = "Only two confirmed peat contacts; each row has its own shared colour scale.",
    caption = "Black dots: confirmed contacts. Depth and carbon are updated separately; the carbon map is not the local curve extrapolated.\nLi depth/carbon priors are dependent. These maps do not provide a validated joint peat-carbon inventory.",
    theme = theme(plot.title = element_text(size = 19, face = "bold"),
      plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 11, hjust = 0)))
ggsave(file.path(out, "li_map_comparison.png"), li, width = 11, height = 12, dpi = 160, bg = "white")
ggsave(file.path(out, "li_map_comparison.pdf"), li, width = 11, height = 12, device = cairo_pdf)

write.csv(data.frame(method = c("Sothe prior", "Spatial core update", "Area-based update", "IDW update"),
  n_cells = 3873, map_mean_kg_C_m2 = c(31.38781654733153, 29.975413725505344, 13.17771554790741, 13.987392140023475)),
  file.path(out, "sothe_summary.csv"), row.names = FALSE)
