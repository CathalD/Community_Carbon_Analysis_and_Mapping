# Translate raw measured FS layers to the workshop's field workbook schema.
library(dplyr)
layers <- read.csv("data/field/core_layers.csv", check.names = FALSE)
locations <- read.csv("data/field/core_locations.csv", check.names = FALSE)
stopifnot(!anyDuplicated(locations$plot_id), setequal(layers$plot_id, locations$plot_id),
  all(is.finite(layers$soc)), all(layers$soc >= 0 & layers$soc <= 1000),
  all(is.finite(layers$bulk_density) & layers$bulk_density > 0),
  all(layers$coarse_frag >= 0 & layers$coarse_frag <= 1),
  all(layers$depth_to > layers$depth_from & layers$depth_from >= 0))
layers <- layers[order(layers$plot_id, layers$depth_from), ]
core_sheet <- data.frame(`Core ID` = locations$plot_id, `Plot ID` = locations$plot_id,
  Latitude = locations$latitude, Longitude = locations$longitude, check.names = FALSE)
# Raw SOC is g/kg; workshop SOC is percent. Stock is kg C/m2.
sample_sheet <- data.frame(`Plot ID` = layers$plot_id, `Core ID` = layers$plot_id,
  `Top depth (cm)` = layers$depth_from, `Bottom depth (cm)` = layers$depth_to,
  `Bulk density (g/cm³)` = layers$bulk_density, `Organic carbon (%)` = layers$soc / 10,
  `Stock, full thickness (kg C/m²)` = layers$soc * layers$bulk_density *
    (1 - layers$coarse_frag) * (layers$depth_to - layers$depth_from) / 100,
  check.names = FALSE)
sample_sheet[["Cumulative carbon to base (kg C/m²)"]] <- ave(
  sample_sheet[["Stock, full thickness (kg C/m²)"]], sample_sheet[["Core ID"]], FUN = cumsum)
