
#==========================================

# ----- Part 4 ------------------------

#======================================


# ===============
# --- Step 1 ----
#================
#Import all the maps that were created so far

cores <- read_csv(
  "Draft1/data/mapping_data.csv",
  show_col_types = FALSE
)

# Comparing map outputs and prior data
library(sf)
st_layers("Draft1/outputs/area_based.gpkg")
st_layers("Draft1/outputs/aoi_polygon.gpkg")
# Step 1 - Load current maps:
aoi <- vect(st_read("Draft1/outputs/aoi_polygon.gpkg"))
# Area-based extrapolation
ab_0_30 <- st_read(
  area_based_gpkg,
  layer = "area_based_0_30"
)
ab_0_100 <- st_read(
  area_based_gpkg,
  layer = "area_based_0_100"
)
# IDW
idw_0_30 <- rast("Draft1/outputs/IDW_0_30.tif")
idw_0_100 <- rast("Draft1/outputs/IDW_0_100.tif")

# Raster-based modelling
# Leave alone for now


# ===============
# --- Step 2 ----
#================
# Area-based vs IDW
# A) table for mean + uncertainty

# Summary Statistics for area-based (The actual field data is used here)
mean_ab <- mean(cores$`C Stock 0-100cm`)
sd_ab   <- sd(cores$`C Stock 0-100cm`)
n_ab    <- nrow(cores)
se_ab <- sd_ab / sqrt(n_ab)
cv_ab <- 100 * sd_ab / mean_ab

ci_ab <- c(
  lower = mean_ab - 1.96 * se_ab,
  upper = mean_ab + 1.96 * se_ab
)


# Summary Statistcs for IDW

idw_vals <- na.omit(values(idw_0_100))
mean_idw <- mean(idw_vals)
sd_idw <- sd(idw_vals)
cv_idw <- 100 * sd_idw / mean_idw
ci <- function(x){
  m <- mean(x)
  se <- sd(x)/sqrt(length(x))
  c(
    lower = m - 1.96*se,
    upper = m + 1.96*se
  )
}

ci_ab  <- ci(cores$`C Stock 0-100cm`)
ci_idw <- ci(idw_vals)

summary_tbl <- data.frame(
  
  Method = c("Area-based","IDW"),
  
  Mean = c(mean_ab,mean_idw),
  
  SD = c(sd_ab,sd_idw),
  
  CV = c(cv_ab,cv_idw),
  
  Lower95 = c(ci_ab[1],ci_idw[1]),
  
  Upper95 = c(ci_ab[2],ci_idw[2])
  
)

summary_tbl



# Comparing AB + IDW:
# Harmonize data
ab_raster <- rasterize(
  vect(ab_0_100),
  idw_0_100,
  field = "mean_stock_kg_m2"
)

ab_vals  <- values(ab_raster)
idw_vals <- values(idw_0_100)
compare_df <- data.frame(
  area_based = ab_vals,
  idw = idw_vals
)
compare_df <- na.omit(compare_df)

ab_vals  <- as.numeric(values(ab_raster)[,1])
idw_vals <- as.numeric(values(idw_0_100)[,1])

compare_df <- data.frame(
  area_based = ab_vals,
  idw = idw_vals
)

compare_df <- na.omit(compare_df)
cor <- cor(compare_df$area_based, compare_df$idw, method = "spearman")
cor
rmse <- sqrt(
  mean(
    (compare_df$idw -
       compare_df$area_based)^2
  )
)

rmse



# B) map of Area-based - IDW

diff_raster <- idw_0_100 - ab_raster

global(diff_raster, "mean", na.rm = TRUE)
global(diff_raster, "sd", na.rm = TRUE)
global(diff_raster, "min", na.rm = TRUE)
global(diff_raster, "max", na.rm = TRUE)


plot(
  diff_raster,
  col = hcl.colors(
    25,
    "Blue-Red 3",
    rev = TRUE
  ),
  main = "IDW - Area-Based"
)

hist(
  values(diff_raster),
  breaks = 30,
  main = "Distribution of Differences",
  xlab = "IDW - Area-Based"
)


# C)Map selection (which is most parsimonious)


# Prior data:
# 1 - Sothe et al. (2022) - 250m
# Soil Carbon stocks to 30 cm
# Soil Carbon to 100 cm + uncertianty

prior_MM_0_100 <- rast("Y:/GIS_CORE/CORE/Canada/Physical/Carbon/McMasterCarbonData/Version5.0/McMaster_WWFCanada_soil_carbon1m_250m_kg-m2_version3.0.tif")
prior_MM_sd <- rast("Y:/GIS_CORE/CORE/Canada/Physical/Carbon/McMasterCarbonData/Version5.0/McMaster_WWFCanada_soil_carbon1m_uncertainty_250m_kg-m2_version3.0.tif")

# I need rgee for the next two
# 2 - WOSIS - SoilsGrid database - 250 m
# Soil Carbon stock present in 5 image bands in GEE
# 0-5, 5 - 15, 15 - 30, 30- 60, 60 - 100 cm
# 3 - Agriculture Canada soil map - 100m
# 5 bands for different depths

# For each, list out average carbon stock across AOI

# Harmonize Data
library(terra)
#  Clip to AOI
crs(prior_MM_0_100)
st_crs(aoi)
aoi <- project(aoi, crs(prior_MM_0_100))

r_mm_clip <- prior_MM_0_100 |>
  crop(aoi) |>
  mask(aoi)

r_sd_clip <- prior_MM_sd |>
  crop(aoi) |>
  mask(aoi)
plot(r_mm_clip)

# Mathc crs
r_mm_proj <- project(
  r_mm_clip,
  idw_0_100,
  method = "bilinear"
)

r_sd_proj <- project(
  r_sd_clip,
  idw_0_100,
  method = "bilinear"
)


core_prior <- terra::extract(
  r_mm_proj,
  vect(core_proj)
)

core_prior_df <- data.frame(
  observed = cores$`C Stock 0-100cm`,
  prior = core_prior[,2]
)

core_prior_df
summary(core_prior_df)

mean_obs <- mean(
  core_prior_df$observed,
  na.rm = TRUE
)

mean_prior <- mean(
  core_prior_df$prior,
  na.rm = TRUE
)

sd_obs <- sd(
  core_prior_df$observed,
  na.rm = TRUE
)

sd_prior <- sd(
  core_prior_df$prior,
  na.rm = TRUE
)

core_summary <- data.frame(
  Dataset = c("Observed Cores", "McMaster Prior"),
  Mean = c(mean_obs, mean_prior),
  SD = c(sd_obs, sd_prior)
)

core_summary


# Difference Between Area-based and Prior
diff_ab_prior <- ab_raster - r_mm_proj
compareGeom(
  ab_raster,
  r_mm_proj,
  stopOnError = FALSE
)

global(
  diff_ab_prior,
  "mean",
  na.rm = TRUE
)

pct_ab_prior <- (
  (ab_raster - r_mm_proj) /
    r_mm_proj
) * 100

plot(
  pct_ab_prior,
  main = "Area-Based vs Prior (%)"
)

vals <- values(pct_ab_prior)

agree20 <- 100 *
  sum(abs(vals) <= 20, na.rm = TRUE) /
  sum(!is.na(vals))

agree20


ab_mean <- global(
  ab_raster,
  "mean",
  na.rm = TRUE
)[1,1]

prior_mean <- global(
  r_mm_proj,
  "mean",
  na.rm = TRUE
)[1,1]

ab_sd <- sd(
  values(ab_raster),
  na.rm = TRUE
)

prior_sd <- sd(
  values(r_mm_proj),
  na.rm = TRUE
)

data.frame(
  Dataset = c("Area-Based", "McMaster Prior"),
  Mean = c(ab_mean, prior_mean),
  SD = c(ab_sd, prior_sd)
)




# Repeat for IDW - Prior Comparison






# ===================================
# ----- Exports -------------------
# ================================

#Export priors to output folder for next step

writeRaster(
  r_mm_proj,
  "Draft1/outputs/prior_MM_0_100.tif",
  overwrite = TRUE
)

writeRaster(
  r_sd_proj,
  "Draft1/outputs/prior_MM_sd.tif",
  overwrite = TRUE
)





