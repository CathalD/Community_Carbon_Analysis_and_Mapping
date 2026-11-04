#===============================================================

# -------- Part 2 - Begin Spatial Mapping + Analysis

#===============================================================

#### notes - Config crs so I dont have to add manually###

#If packages not installed - please do so one time:

#install.packages("sf")
#install.packages("terra")
#install.packages("gstat")
#install.packages("patchwork")

# ===================

# ---- Step 1 ------

# Read in the field data from Part 1:

library(sf)
library(readr)
library(dplyr)
library(ggplot2)
#create a place to store outputs
dir.create("outputs/script2", recursive = TRUE, showWarnings = FALSE)

#----------------
# Read input files
#-------------

cores <- read_csv(
  "outputs/script1/mapping_data.csv",
  show_col_types = FALSE
)

aoi <- st_read(
  "data/boundary/aoi.geojson",
  quiet = TRUE
)
# View imports
glimpse(cores)
head(cores)
st_geometry_type(aoi)
st_crs(aoi)
# Convert cores to spatial points
cores <- cores %>%
  rename(
    core_id = `Core ID`,
    latitude = Latitude,
    longitude = Longitude,
    stock_0_30 = `C Stock 0-30cm`,
    stock_0_50 = `C Stock 0-50cm`,
    stock_0_100 = `C Stock 0-100cm`,
    stock_Total = `Total Carbon Stock`
  )
core_sf <- st_as_sf(
  cores,
  coords = c("longitude", "latitude"),
  crs = 4326
)
print(core_sf)
# GIS housekeeping - Get CRS to same projection
aoi_proj <- st_transform(aoi, 3978)
core_proj <- st_transform(core_sf, 3978)

st_write(aoi_proj, "outputs/script2/aoi_polygon.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
# Inspect the AOI and Points
ggplot() +
  geom_sf(data = aoi_proj,
          fill = "grey95",
          color = "black") +
  geom_sf(data = core_proj,
          aes(color = stock_0_30),
          size = 3) +
  scale_color_viridis_c() +
  theme_minimal()

# ========================
# ---- Step 2A ----------
# Simple Area-based extrapolation

library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)

# -----------------------------------
# AOI area calculations/ Zonal Stats

aoi_area_m2 <- as.numeric(st_area(aoi_proj))
aoi_area_ha <- aoi_area_m2 / 10000
aoi_area_km2 <- aoi_area_m2 / 1e6

cat("AOI Area:", round(aoi_area_km2,2), "km²\n")

# ------------------------
# Calculate mean
mean_30 <- mean(core_sf$stock_0_30, na.rm = TRUE)
mean_50 <- mean(core_sf$stock_0_50, na.rm = TRUE)
mean_100 <- mean(core_sf$stock_0_100, na.rm = TRUE)
mean_total <- mean(core_sf$stock_Total, na.rm = TRUE)

# -----------------------------------
# Calculate AOI totals
total_30_tC    <- mean_30    * aoi_area_m2 / 1000
total_50_tC    <- mean_50    * aoi_area_m2 / 1000
total_100_tC   <- mean_100   * aoi_area_m2 / 1000
total_total_tC <- mean_total * aoi_area_m2 / 1000

# ---------------------------
# Summary table
area_summary <- data.frame(
  Depth = c(
    "0-30 cm",
    "0-50 cm",
    "0-100 cm",
    "Full Profile"
  ),
  Mean_Stock_kg_m2 = c(
    mean_30,
    mean_50,
    mean_100,
    mean_total
  ),
  Total_Carbon_tC = c(
    total_30_tC,
    total_50_tC,
    total_100_tC,
    total_total_tC
  )
)

print(area_summary)

stock_vars <- c(
  "stock_0_30",
  "stock_0_50",
  "stock_0_100",
  "stock_Total"
)

titles <- c(
  "0-30 cm",
  "0-50 cm",
  "0-100 cm",
  "Full Profile"
)

# Means corresponding to above variables

mean_values <- c(
  mean_30,
  mean_50,
  mean_100,
  mean_total
)

#Colour

all_values <- c(
  core_sf$stock_0_30,
  core_sf$stock_0_50,
  core_sf$stock_0_100,
  core_sf$stock_Total,
  mean_values
)

common_limits <- range(all_values, na.rm = TRUE)

# -----------------------------------
# Function to build Maps for each depth

create_area_map <- function( # AI helped write this function
    stock_var,
    title,
    mean_value
){
  
  aoi_plot <- aoi_proj
  aoi_plot$mean_stock <- mean_value
  
  ggplot() +
    
    # AOI coloured by mean value
    geom_sf(
      data = aoi_plot,
      aes(fill = mean_stock),
      colour = "black",
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    
    # Core observations
    geom_sf(
      data = core_proj,
      aes(colour = .data[[stock_var]]),
      size = 4
    ) +
    
    scale_fill_viridis_c(
      limits = common_limits,
      option = "C"
    ) +
    
    scale_colour_viridis_c(
      limits = common_limits,
      option = "C",
      name = "Carbon Stock\n(kg C/m²)"
    ) +
    
    labs(
      title = title,
      subtitle = paste0(
        "AOI Mean = ",
        round(mean_value,1),
        " kg C/m²"
      )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# -----------------------------------
# Make 4 maps for each depth and export to folder
plot_list <- Map(
  create_area_map,
  stock_var = stock_vars,
  title = titles,
  mean_value = mean_values
)
plot_list # Viiw atll 4 plots

output_dir <- "outputs/script2"
library(sf)

st_write(
  aoi_proj,
  "outputs/script2/area_based.gpkg",
  delete_dsn = TRUE
)

depth_names <- c(
  "0_30",
  "0_50",
  "0_100",
  "Total"
)

area_based_gpkg <- "outputs/script2/area_based.gpkg"

for(i in seq_along(depth_names)) { #AI was used to help write this function
  
  area_poly <- aoi_proj
  area_poly$depth <- depth_names[i]
  area_poly$method <- "area_based"
  area_poly$mean_stock_kg_m2 <- mean_values[i]
  
  st_write(
    obj = area_poly,
    dsn = area_based_gpkg,
    layer = paste0("area_based_", depth_names[i]),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
}

# Option to view:
st_layers(area_based_gpkg)
# -----------------------------------
# Combine into 2x2 figure
# -----------------------------------

combined_plot <-
  wrap_plots(plot_list, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

combined_plot


# ---------------------------
# Here we build on the area-averaged mean stock,
# adding in a distance weighted function using the 
# IDW function in the package gstat

# The code was directed adapted from this source:

# For documentation see here: https://www.paulamoraga.com/book-spatial/spatial-interpolation-methods.html

# Package: - gstat
# Function: gstat::idw
# Parameters: Decay function + distance weighting

library(sf)
library(terra)

# First build a grid using sf function
# You can define the cell size here
##### NOTE add cellsize as a config #######
grid <- st_make_grid(
  aoi_proj,
  cellsize = 1000,
  square = TRUE
  
)
grid_sf <- st_sf(
  grid_id = seq_along(grid),
  geometry = grid
)
grid_sf <- grid_sf[
  lengths(st_intersects(grid_sf, aoi_proj)) > 0,
]
nrow(grid_sf)
ggplot() +
  geom_sf(data = grid_sf,
          fill = NA,
          colour = "grey80",
          linewidth = 0.2) +
  geom_sf(data = aoi_proj,
          fill = NA,
          colour = "black") +
  geom_sf(data = core_proj,
          colour = "red",
          size = 3)

# Assign a point to middle of each grid
grid_points <- st_centroid(grid_sf)
coords <- st_coordinates(grid_points)
grid_points$x <- coords[,1]
grid_points$y <- coords[,2]
head(grid_points)
head(core_proj)
count(grid_points)
ggplot() +
  geom_sf(data = grid_points,
          size = 0.0001,
          colour = "blue") +
  geom_sf(data = aoi_proj,
          fill = NA,
          colour = "black") +
  geom_sf(data = core_proj,
          colour = "red",
          size = 3)

# Model 0 - 30 cm
library(gstat)

pred_points = grid_points

# Create polygons for grouping AOI by closest core
library(terra)
vor <- terra::voronoi(
  x = terra::vect(core_proj),
  bnd = terra::vect(aoi_proj)
)
vor_sf <- st_as_sf(vor)
ggplot() +
  geom_sf(
    data = vor_sf,
    aes(fill = stock_0_30)
  ) +
  geom_sf(
    data = core_proj,
    colour = "black",
    size = 3
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

closest_pred <- st_intersection(
  pred_points,
  vor_sf
)
head(closest_pred$stock_0_30)


###### Here resolution is also defined and must be same as CELLSIZE #####

r <- rast(
  vect(aoi_proj),
  resolution = 1000
)

closest_raster <- terra::rasterize(
  vect(closest_pred),
  r,
  field = "stock_0_30"
)

closest_raster <- terra::rasterize(
  vect(closest_pred),
  r,
  field = "stock_0_30"
)
plot(closest_raster)   


# Lets run the IDW prediction for stocks 0-30cm first
summary(closest_pred$stock_0_30)

IDW_30 <- gstat(formula = stock_0_30 ~ 1,
                locations = core_proj,
                nmax = nrow(core_proj), # each core has 5 neighbours
                set = list(idp = 1)) # This controls the strength of the distance as an exponetial function, setting this to 0 = area-based extrapoloation

IDW_30_pred <- predict(
  IDW_30,
  pred_points
)

head(IDW_30_pred)

summary(core_proj$stock_0_30)
summary(IDW_30_pred$var1.pred)
ggplot() +
  geom_sf(
    data = IDW_30_pred,
    aes(colour = var1.pred)
  ) +
  scale_colour_viridis_c() +
  theme_minimal()


#------------------------

#Rasterize the outputs:

IDW_30_pred$pred<- IDW_30_pred$var1.pred

IDW_30_raster <- terra::rasterize(
  vect(IDW_30_pred),
  r,
  field = "pred"
)
plot(IDW_30_raster)


# Create a carbon map and export
idw_df <- as.data.frame(
  IDW_30_raster,
  xy = TRUE,
  na.rm = TRUE
)

names(idw_df)[3] <- "stock"

ggplot() +
  
  geom_raster(
    data = idw_df,
    aes(
      x = x,
      y = y,
      fill = stock
    )
  ) +
  
  geom_sf(
    data = aoi_proj,
    fill = NA,
    colour = "black"
  ) +
  
  geom_sf(
    data = core_proj,
    colour = "white",
    fill = "black",
    shape = 21,
    size = 3
  ) +
  
  scale_fill_viridis_c(
    name = "SOC Stock\n(kg C/m²)"
  ) +
  
  coord_sf() +
  
  theme_minimal() +
  
  labs(
    title = "0-30 cm Carbon Stock",
    subtitle = "Inverse Distance Weighting (β = 1)"
  )

# Now repeat for the other depths and stack as bands for same geotiff for export


# Here we are going to run IDW methods for the 4 differnet depth intervals
# WE are going to measure the variance from each test using Leave-on-out core validatio
# and record the RMSE -

# For each depth, the best performing parameter settings will be saved and exported for
# Further Analysis

# Parameters we modify are: 
#       1. the distance power, represented by idp in idw function
#       2. The number fo neighbors (nMax)
#       3. (not tested yet) Setting max distance (area-of-applicability)

library(gstat)

idw_loocv <- function(data, variable, idp, nmax) { # AI helped write this function
  predicted <- numeric(nrow(data))
  for (i in seq_len(nrow(data))) {
    training <- data[-i, ]
    model <- gstat::gstat(
      formula = stats::as.formula(paste(variable, "~ 1")),
      locations = training,
      nmax = nmax,
      set = list(idp = idp)
    )
    predicted[i] <- stats::predict(
      model, data[i, ], debug.level = 0
    )$var1.pred
  }
  sqrt(mean((data[[variable]] - predicted)^2))
}

tune_depth <- function(variable, depth_label) {
  data <- core_proj[is.finite(core_proj[[variable]]), ]
  parameters <- expand.grid(
    idp = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5), # NOTE - When re-running - Tune this loop to what fits best with the data you are using
    nmax = c(2, 3, 4, 5)
  )
  # A withheld core must never count as a training neighbour.
  parameters <- parameters[parameters$nmax <= nrow(data) - 1L, ]
  parameters$RMSE <- vapply(seq_len(nrow(parameters)), function(i) { # AI helped write this function
    idw_loocv(data, variable, parameters$idp[i], parameters$nmax[i])
  }, numeric(1))
  
  best <- parameters[which.min(parameters$RMSE), , drop = FALSE]
  model <- gstat::gstat(
    formula = stats::as.formula(paste(variable, "~ 1")),
    locations = data,
    nmax = best$nmax,
    set = list(idp = best$idp)
  )
  prediction <- stats::predict(model, pred_points, debug.level = 0)
  prediction$pred <- prediction$var1.pred
  
  parameters$Depth <- depth_label
  parameters$Variable <- variable
  best$Depth <- depth_label
  best$Variable <- variable
  best$N_Cores <- nrow(data)
  list(best = best, parameters = parameters, prediction = prediction)
}

depth_variables <- c(
  "0-30 cm" = "stock_0_30",
  "0-50 cm" = "stock_0_50",
  "0-100 cm" = "stock_0_100",
  "Full Profile" = "stock_Total"
)

depth_results <- lapply(names(depth_variables), function(label) {
  tune_depth(depth_variables[[label]], label)
})
names(depth_results) <- unname(depth_variables)

# Preserve the object names expected by the rest of Script_2.R.
pred_30 <- depth_results[["stock_0_30"]]$prediction
pred_50 <- depth_results[["stock_0_50"]]$prediction
pred_100 <- depth_results[["stock_0_100"]]$prediction
pred_total <- depth_results[["stock_Total"]]$prediction

best_by_depth <- do.call(rbind, lapply(depth_results, function(x) x$best))
best_by_depth <- best_by_depth[, c("Depth", "Variable", "N_Cores", "idp", "nmax", "RMSE")]
rownames(best_by_depth) <- NULL
tuning_by_depth <- do.call(rbind, lapply(depth_results, function(x) x$parameters))
rownames(tuning_by_depth) <- NULL
print(best_by_depth)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(best_by_depth, "outputs/script2/idw_best_parameters_by_depth.csv", row.names = FALSE)
write.csv(tuning_by_depth, "outputs/script2/idw_tuning_by_depth.csv", row.names = FALSE)

# Optional: compare the tuning surfaces for all four depths.
optimization_plot <- ggplot2::ggplot( # AI helped write this visualization function
  tuning_by_depth,
  ggplot2::aes(x = factor(nmax), y = factor(idp), fill = RMSE)
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = round(RMSE, 2)), size = 3) +
  ggplot2::scale_fill_viridis_c(direction = -1) +
  ggplot2::facet_wrap(~ Depth) +
  ggplot2::labs(x = "Neighbours (NMAX)", y = "Distance power (IDP)",
                fill = "Tuning RMSE\n(kg C/m²)") +
  ggplot2::theme_minimal()
print(optimization_plot)




# -------------
library(ggplot2)

ggplot(
  tuning_by_depth,
  aes(
    x = idp,
    y = RMSE,
    colour = factor(nmax),
    group = nmax
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Depth, scales = "free_y") +
  labs(
    x = "IDW Power (IDP)",
    y = "LOOCV RMSE",
    colour = "NMAX"
  ) +
  theme_minimal()



# -------------------------------------
# Rasterize predictions


r30 <- rasterize(
  vect(pred_30),
  r,
  field = "pred"
)

r50 <- rasterize(
  vect(pred_50),
  r,
  field = "pred"
)

r100 <- rasterize(
  vect(pred_100),
  r,
  field = "pred"
)

rtotal <- rasterize(
  vect(pred_total),
  r,
  field = "pred"
)

# -------------------------------------
# Name layers


names(r30) <- "Stock_0_30cm"
names(r50) <- "Stock_0_50cm"
names(r100) <- "Stock_0_100cm"
names(rtotal) <- "Stock_Total"

# -------------------------------------
# Create raster of all layers and export

carbon_stack <- c(
  r30,
  r50,
  r100,
  rtotal
)

writeRaster(
  r30,
  "outputs/script2/idw_0_30cm.tif",
  overwrite = TRUE
)

writeRaster(
  r50,
  "outputs/script2/idw_0_50cm.tif",
  overwrite = TRUE
)

writeRaster(
  r100,
  "outputs/script2/idw_0_100cm.tif",
  overwrite = TRUE
)

writeRaster(
  rtotal,
  "outputs/script2/idw_total.tif",
  overwrite = TRUE
)

# =====================================
# Create maps

raster_to_df <- function(rast_obj){
  
  df <- as.data.frame(
    rast_obj,
    xy = TRUE,
    na.rm = TRUE
  )
  
  names(df)[3] <- "stock"
  
  return(df)
  
}

df30 <- raster_to_df(r30)
df50 <- raster_to_df(r50)
df100 <- raster_to_df(r100)
dfTotal <- raster_to_df(rtotal)

# Common colour scale
all_values <- c(
  values(r30),
  values(r50),
  values(r100),
  values(rtotal)
)

common_limits <- range(
  all_values,
  na.rm = TRUE
)

# -------------------------------------
# Mapping function

create_idw_map <- function(
    raster_df,
    title
){
  
  ggplot() +
    
    geom_raster(
      data = raster_df,
      aes(
        x = x,
        y = y,
        fill = stock
      )
    ) +
    
    geom_sf(
      data = aoi_proj,
      fill = NA,
      colour = "black",
      linewidth = 0.4
    ) +
    
    geom_sf(
      data = core_proj,
      colour = "white",
      fill = "black",
      shape = 21,
      size = 3
    ) +
    
    scale_fill_viridis_c(
      limits = common_limits,
      name = "kg C/m²"
    ) +
    
    coord_sf() +
    
    theme_minimal() +
    
    labs(
      title = title
    )
  
}

# -------------------------------------
# Build maps

p30_idw <- create_idw_map(
  df30,
  "IDW Carbon Stock (0-30 cm)"
)

p50_idw <- create_idw_map(
  df50,
  "IDW Carbon Stock (0-50 cm)"
)

p100_idw <- create_idw_map(
  df100,
  "IDW Carbon Stock (0-100 cm)"
)

pTotal_idw <- create_idw_map(
  dfTotal,
  "IDW Carbon Stock (Full Profile)"
)

# -------------------------------------
# Combine into 1 figure


library(patchwork)

idw_panel <-
  (p30_idw | p50_idw) /
  (p100_idw | pTotal_idw)

idw_panel

# -------------------------------------
# Save figure

ggsave(
  "outputs/script2/idw_maps.png",
  idw_panel,
  width = 12,
  height = 10,
  dpi = 300
)


# =====================================
# Calculate AOI carbon totals

# Attach predictions back to grid polygons

grid_30 <- grid_sf
grid_30$stock <- pred_30$pred

grid_50 <- grid_sf
grid_50$stock <- pred_50$pred

grid_100 <- grid_sf
grid_100$stock <- pred_100$pred

grid_total <- grid_sf
grid_total$stock <- pred_total$pred

# -------------------------------------
# Clip to AOI
# Gives exact boundary areas
# -------------------------------------

grid_30 <- st_intersection(
  grid_30,
  aoi_proj
)

grid_50 <- st_intersection(
  grid_50,
  aoi_proj
)

grid_100 <- st_intersection(
  grid_100,
  aoi_proj
)

grid_total <- st_intersection(
  grid_total,
  aoi_proj
)

# -------------------------------------
# Calculate cell areas
# -------------------------------------

grid_30$area_m2 <- as.numeric(
  st_area(grid_30)
)

grid_50$area_m2 <- as.numeric(
  st_area(grid_50)
)

grid_100$area_m2 <- as.numeric(
  st_area(grid_100)
)

grid_total$area_m2 <- as.numeric(
  st_area(grid_total)
)

# -------------------------------------
# Calculate IDW carbon totals
# -------------------------------------

total_30_idw <- sum(
  grid_30$stock *
    grid_30$area_m2,
  na.rm = TRUE
) / 1000

total_50_idw <- sum(
  grid_50$stock *
    grid_50$area_m2,
  na.rm = TRUE
) / 1000

total_100_idw <- sum(
  grid_100$stock *
    grid_100$area_m2,
  na.rm = TRUE
) / 1000

total_profile_idw <- sum(
  grid_total$stock *
    grid_total$area_m2,
  na.rm = TRUE
) / 1000

# -------------------------------------
# IDW summary
# -------------------------------------

idw_summary <- data.frame(
  Depth = c(
    "0-30 cm",
    "0-50 cm",
    "0-100 cm",
    "Full Profile"
  ),
  IDW_Total_tC = c(
    total_30_idw,
    total_50_idw,
    total_100_idw,
    total_profile_idw
  )
)

print(idw_summary)

# -------------------------------------
# Compare with area-based extrapolation
# -------------------------------------

comparison <- data.frame(
  Depth = c(
    "0-30 cm",
    "0-50 cm",
    "0-100 cm",
    "Full Profile"
  ),
  
  Area_Average_tC = c(
    total_30_tC,
    total_50_tC,
    total_100_tC,
    total_total_tC
  ),
  
  IDW_tC = c(
    total_30_idw,
    total_50_idw,
    total_100_idw,
    total_profile_idw
  )
)

comparison$Difference_tC <-
  comparison$IDW_tC -
  comparison$Area_Average_tC

comparison$Percent_Difference <-
  100 * (
    comparison$IDW_tC -
      comparison$Area_Average_tC
  ) /
  comparison$Area_Average_tC

print(comparison)

# -------------------------------------
# Export summaries
# -------------------------------------

write.csv(
  idw_summary,
  "outputs/script2/idw_carbon_totals.csv",
  row.names = FALSE
)

write.csv(
  comparison,
  "outputs/script2/idw_vs_area_average.csv",
  row.names = FALSE
)

# -------------------------------------
# Quick visualization
# -------------------------------------

library(ggplot2)

ggplot(
  comparison,
  aes(
    x = Depth,
    y = Area_Average_tC
  )
) +
  geom_col(
    fill = "grey70",
    width = 0.7
  ) +
  geom_point(
    aes(y = IDW_tC),
    colour = "red",
    size = 4
  ) +
  labs(
    title = "Carbon Stock Comparison",
    subtitle = "Grey bars = Area Average, Red points = Optimized IDW",
    y = "Total Carbon Stock (t C)",
    x = NULL
  ) +
  theme_minimal()



# ================================
# Here I should bring in Ordinary kriging
# =====================================
# Ordinary Kriging
# =====================================

library(gstat)

# -----------------------
# Fit empirical variogram

#First log transofrm the data

stock_0_30_log = log(core_proj$stock_0_30)

ok_variogram <- variogram(
  stock_0_30_log ~ 1,
  locations = core_proj
)

plot(ok_variogram)

# -----------------------
# Fit theoretical model
# Try Exp, Sph, Gau

ok_model <- fit.variogram(
  ok_variogram,
  model = vgm("Sph")
)

plot(
  ok_variogram,
  model = ok_model
)

print(ok_model)

# -----------------------
# Perform Ordinary Kriging

OK_30 <- krige(
  formula = stock_0_30_log ~ 1,
  locations = core_proj,
  newdata = pred_points,
  model = ok_model
)

head(OK_30)

OK_30$pred <- OK_30$var1.pred

OK_30_raster <- rasterize(
  vect(OK_30),
  r,
  field = "pred"
)

plot(OK_30_raster)

OK_30$variance <- OK_30$var1.var

OK_30_var_raster <- rasterize(
  vect(OK_30),
  r,
  field = "variance"
)

plot(OK_30_var_raster)






