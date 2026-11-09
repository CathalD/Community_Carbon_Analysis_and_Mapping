"""Pasteable Colab cells for Canada-wide external-core covariate sampling.

This section does not resubmit the completed AOI raster exports. In a fresh
runtime, run the notebook's import, authentication/configuration, and helper
cells first (Sections 0 and 2), then paste/run the sections below.

The input is the Canada-only table written by
scripts/external_profiles/06_build_external_core_table.R. Output is split by
broad ecosystem group and, when necessary, into bounded batches. The union of
the output CSV files is the complete Canada sample selected below.
"""


# %% --------------------------------------------------------------------------
# 9. LOAD THE CANADA-ONLY CORE TABLE
# -----------------------------------------------------------------------------

import math
import os
import re

import ee
import pandas as pd

# Reuse cached Colab credentials when available. This initializes Earth Engine
# but does not start either the raster or point export tasks.
try:
    ee.Initialize(project="northstarlabs")
except Exception:
    ee.Authenticate()
    ee.Initialize(project="northstarlabs")

required_notebook_objects = [
    "CONFIG",
    "FOCAL_PROJECTION",
    "fixed_calendar_window",
    "mask_s2_clouds",
    "add_spring_wetness_indices",
    "add_summer_greenness_indices",
    "add_focal_stats",
]
missing_notebook_objects = [
    name for name in required_notebook_objects if name not in globals()
]
if missing_notebook_objects:
    raise RuntimeError(
        "Run the notebook import, configuration and helper cells first. "
        "Do not run the raster-export cell. Missing definitions: "
        + ", ".join(missing_notebook_objects)
    )

try:
    from google.colab import files
    IN_COLAB = True
except ImportError:
    IN_COLAB = False

CORE_TABLE_LOCAL_PATH = "/content/external_carbon_cores_canada_for_gee.csv"

if not os.path.exists(CORE_TABLE_LOCAL_PATH):
    if IN_COLAB:
        uploaded = files.upload()
        CORE_TABLE_LOCAL_PATH = next(iter(uploaded))
    else:
        raise FileNotFoundError(
            f"{CORE_TABLE_LOCAL_PATH} not found. Upload the Canada-only CSV "
            "or change CORE_TABLE_LOCAL_PATH before running this cell."
        )

required_core_cols = {
    "core_uid",
    "database_code",
    "latitude",
    "longitude",
    "country",
    "source_ecosystem",
    "broad_ecosystem_group",
    "has_any_complete_interval",
}
cores_df = pd.read_csv(CORE_TABLE_LOCAL_PATH, low_memory=False)
missing_core_cols = required_core_cols - set(cores_df.columns)
if missing_core_cols:
    raise ValueError(
        "Rebuild the core tables before continuing. Missing columns: "
        f"{sorted(missing_core_cols)}"
    )

cores_df = cores_df[
    cores_df["country"].astype(str).str.casefold().eq("canada")
].copy()
cores_df = cores_df.dropna(subset=["latitude", "longitude"])
cores_df = cores_df[
    cores_df["latitude"].between(-90, 90)
    & cores_df["longitude"].between(-180, 180)
].drop_duplicates("core_uid")

print(f"Loaded {len(cores_df):,} unique Canadian cores")
print(cores_df["database_code"].value_counts().sort_index())


# %% --------------------------------------------------------------------------
# 10. SELECT CARBON-ELIGIBLE CORES AND PLAN ECOSYSTEM EXPORTS
# -----------------------------------------------------------------------------

# Default pilot selection: Riley and CanPeat. Set this to None to retain all
# database sources in the Canada-only input table.
DATABASE_FILTER = ["RILEY", "CANPEAT"]

# Keep every selected core in the extraction table and use its carbon coverage
# fields during modelling. Set True when only response-complete cores are wanted.
REQUIRE_COMPLETE_CARBON_INTERVAL = False

# Leave as None to export every broad ecosystem group. To run only selected
# groups, supply exact names printed by the export plan below.
ECOSYSTEM_GROUP_FILTER = None
# Example:
# ECOSYSTEM_GROUP_FILTER = ["peatland_bog", "peatland_fen"]

# Optional spatial restriction. False keeps all selected Canadian cores. True
# keeps only cores inside the AOI bounding box. Increase the padding only when
# you intentionally want a wider rectangle around the uploaded AOI.
AOI_bounds = False
AOI_BOUNDS_PAD_DEG = 0.0

selected_cores_df = cores_df.copy()
if DATABASE_FILTER is not None:
    selected_cores_df = selected_cores_df[
        selected_cores_df["database_code"].isin(DATABASE_FILTER)
    ]
if REQUIRE_COMPLETE_CARBON_INTERVAL:
    selected_cores_df = selected_cores_df[
        selected_cores_df["has_any_complete_interval"].fillna(False)
    ]
if ECOSYSTEM_GROUP_FILTER is not None:
    selected_cores_df = selected_cores_df[
        selected_cores_df["broad_ecosystem_group"].isin(ECOSYSTEM_GROUP_FILTER)
    ]

aoi_filter_bounds = None
if AOI_bounds:
    source_aoi = ee.FeatureCollection(CONFIG["aoi_asset"]).geometry(maxError=1)
    source_aoi_coords = (
        source_aoi.bounds(maxError=1).coordinates().get(0).getInfo()
    )
    source_aoi_lons = [point[0] for point in source_aoi_coords]
    source_aoi_lats = [point[1] for point in source_aoi_coords]
    aoi_lon_min = min(source_aoi_lons) - AOI_BOUNDS_PAD_DEG
    aoi_lon_max = max(source_aoi_lons) + AOI_BOUNDS_PAD_DEG
    aoi_lat_min = min(source_aoi_lats) - AOI_BOUNDS_PAD_DEG
    aoi_lat_max = max(source_aoi_lats) + AOI_BOUNDS_PAD_DEG
    selected_cores_df = selected_cores_df[
        selected_cores_df["longitude"].between(aoi_lon_min, aoi_lon_max)
        & selected_cores_df["latitude"].between(aoi_lat_min, aoi_lat_max)
    ]
    aoi_filter_bounds = ee.Geometry.Rectangle(
        [aoi_lon_min, aoi_lat_min, aoi_lon_max, aoi_lat_max],
        proj="EPSG:4326",
        geodesic=False,
    )
    print(
        "AOI bounding-box filter enabled: "
        f"{aoi_lon_min:.6f}, {aoi_lat_min:.6f}, "
        f"{aoi_lon_max:.6f}, {aoi_lat_max:.6f} "
        f"(padding {AOI_BOUNDS_PAD_DEG} degrees)"
    )
else:
    print("AOI bounding-box filter disabled: retaining the Canada-wide selection")
selected_cores_df = selected_cores_df.reset_index(drop=True)

ecosystem_export_plan = (
    selected_cores_df.groupby(
        ["broad_ecosystem_group", "database_code"], dropna=False
    )
    .size()
    .rename("n_cores")
    .reset_index()
    .sort_values(["broad_ecosystem_group", "database_code"])
)

print(f"Selected {len(selected_cores_df):,} Canadian cores for sampling")
print(ecosystem_export_plan.to_string(index=False))

# A small, separate spreadsheet describing the ecosystem-specific exports.
ECOSYSTEM_PLAN_LOCAL_PATH = "/content/canada_core_ecosystem_export_plan.csv"
ecosystem_export_plan.to_csv(ECOSYSTEM_PLAN_LOCAL_PATH, index=False)
print(f"Wrote ecosystem export plan: {ECOSYSTEM_PLAN_LOCAL_PATH}")


# %% --------------------------------------------------------------------------
# 11. BUILD THE CANADA-WIDE POINT-SAMPLING STACK
# -----------------------------------------------------------------------------
# This rebuilds the same 119 predictors as a lazy Canada-wide computation. It
# does not start or repeat the existing AOI raster exports. Earth Engine only
# evaluates the image at the selected point locations when the table tasks run.

CANADA_BOUNDS = ee.Geometry.Rectangle(
    [-141.1, 41.0, -52.5, 83.2], proj="EPSG:4326", geodesic=False
)
POINT_PROCESSING_BOUNDS = (
    aoi_filter_bounds if AOI_bounds else CANADA_BOUNDS
)

point_year_filter = ee.Filter.calendarRange(
    CONFIG["years"][0], CONFIG["years"][-1], "year"
)
point_spring_filter = fixed_calendar_window(
    *CONFIG["spring_window"][0], *CONFIG["spring_window"][1]
)
point_summer_filter = fixed_calendar_window(
    *CONFIG["summer_window"][0], *CONFIG["summer_window"][1]
)

point_s2_collection = (
    ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(POINT_PROCESSING_BOUNDS)
    .filter(point_year_filter)
    .filter(
        ee.Filter.lt(
            "CLOUDY_PIXEL_PERCENTAGE", CONFIG["s2_cloud_threshold"]
        )
    )
)

point_s2_spring = (
    point_s2_collection.filter(point_spring_filter)
    .map(mask_s2_clouds)
    .map(add_spring_wetness_indices)
)
point_s2_summer = (
    point_s2_collection.filter(point_summer_filter)
    .map(mask_s2_clouds)
    .map(add_summer_greenness_indices)
)

point_spring_count = point_s2_spring.select("NDMI_spring").count().rename(
    "S2_spring_count"
)
point_summer_count = point_s2_summer.select("NDVI_summer").count().rename(
    "S2_summer_count"
)

point_spring_composite = (
    point_s2_spring.select(["NDMI_spring", "TCW_spring"])
    .median()
    .setDefaultProjection(FOCAL_PROJECTION)
    .updateMask(
        point_spring_count.gte(CONFIG["min_s2_observations"])
    )
)
point_summer_composite = (
    point_s2_summer.select(["NDVI_summer", "EVI_summer"])
    .median()
    .setDefaultProjection(FOCAL_PROJECTION)
    .updateMask(
        point_summer_count.gte(CONFIG["min_s2_observations"])
    )
)

point_dw_collection = (
    ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
    .filterBounds(POINT_PROCESSING_BOUNDS)
    .filter(point_year_filter)
    .filter(point_summer_filter)
)
point_dw_count = point_dw_collection.select("trees").count().rename(
    "DW_summer_count"
)
point_dw_probs = (
    point_dw_collection.select(
        ["trees", "shrub_and_scrub", "flooded_vegetation"]
    )
    .mean()
    .rename(["DW_tree_prob", "DW_shrub_prob", "DW_floodedveg_prob"])
    .setDefaultProjection(FOCAL_PROJECTION)
)

point_focal_source_bands = [
    "NDMI_spring",
    "TCW_spring",
    "NDVI_summer",
    "EVI_summer",
    "DW_tree_prob",
    "DW_shrub_prob",
    "DW_floodedveg_prob",
]
point_focal_source = (
    point_spring_composite.addBands(point_summer_composite)
    .addBands(point_dw_probs)
)
point_focal_stats = add_focal_stats(
    point_focal_source,
    point_focal_source_bands,
    CONFIG["focal_radii_m"],
)

point_focal_band_names = []
for radius in CONFIG["focal_radii_m"]:
    point_focal_band_names.extend(
        [
            f"{band}_focalMean_{radius}m"
            for band in point_focal_source_bands
        ]
    )
    point_focal_band_names.extend(
        [
            f"{band}_focalStd_{radius}m"
            for band in point_focal_source_bands
        ]
    )

point_surface_context = point_focal_source.addBands(point_focal_stats).select(
    point_focal_source_bands + point_focal_band_names
)

point_alpha_source_bands = [f"A{axis:02d}" for axis in range(64)]
point_alpha_output_bands = [f"AE_{axis:02d}" for axis in range(64)]
point_alpha_start = f"{CONFIG['alphaearth_year']}-01-01"
point_alpha_end = f"{CONFIG['alphaearth_year'] + 1}-01-01"
point_embeddings = (
    ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")
    .filterBounds(POINT_PROCESSING_BOUNDS)
    .filterDate(point_alpha_start, point_alpha_end)
    .select(point_alpha_source_bands)
    .mosaic()
    .rename(point_alpha_output_bands)
)

point_terrain_band_names = []
point_terrain_stack = None
if CONFIG["include_terrain_context"]:
    point_dem_collection = (
        ee.ImageCollection("NRCan/CDEM")
        .filterBounds(POINT_PROCESSING_BOUNDS)
        .select("elevation")
    )
    point_dem_first = ee.Image(point_dem_collection.first())
    point_dem = (
        point_dem_collection.mosaic()
        .setDefaultProjection(point_dem_first.projection())
        .rename("elevation_m")
    )
    point_slope = ee.Terrain.slope(point_dem).rename("slope_deg")
    point_aspect_rad = ee.Terrain.aspect(point_dem).multiply(
        3.141592653589793 / 180
    )
    point_eastness = point_aspect_rad.sin().rename("eastness")
    point_northness = point_aspect_rad.cos().rename("northness")
    point_terrain_kernel = ee.Kernel.circle(
        radius=CONFIG["terrain_context_radius_m"], units="meters"
    )
    point_local_mean = point_dem.reduceNeighborhood(
        ee.Reducer.mean(), point_terrain_kernel
    )
    point_relative_elevation = point_dem.subtract(point_local_mean).rename(
        f"relative_elevation_{CONFIG['terrain_context_radius_m']}m"
    )
    point_elevation_sd = point_dem.reduceNeighborhood(
        ee.Reducer.stdDev(), point_terrain_kernel
    ).rename(f"elevation_sd_{CONFIG['terrain_context_radius_m']}m")
    point_terrain_stack = ee.Image.cat(
        [
            point_dem,
            point_slope,
            point_eastness,
            point_northness,
            point_relative_elevation,
            point_elevation_sd,
        ]
    )
    point_terrain_band_names = [
        "elevation_m",
        "slope_deg",
        "eastness",
        "northness",
        f"relative_elevation_{CONFIG['terrain_context_radius_m']}m",
        f"elevation_sd_{CONFIG['terrain_context_radius_m']}m",
    ]

point_stack_components = [point_surface_context, point_embeddings]
if point_terrain_stack is not None:
    point_stack_components.append(
        point_terrain_stack.select(point_terrain_band_names)
    )
point_covariate_stack = ee.Image.cat(point_stack_components)

point_qa_stack = ee.Image.cat(
    [
        point_spring_count,
        point_summer_count,
        point_dw_count,
        point_surface_context.mask().reduce(ee.Reducer.min()).rename(
            "surface_context_complete"
        ),
        point_embeddings.mask().reduce(ee.Reducer.min()).rename(
            "embeddings_available"
        ),
    ]
).unmask(0)

# Sentinel values preserve every selected core even if one predictor is masked.
# Use the five QA fields to diagnose coverage before modelling.
MISSING_PREDICTOR_VALUE = -9999
point_sampling_image = point_covariate_stack.unmask(
    MISSING_PREDICTOR_VALUE
).addBands(point_qa_stack)

point_expected_band_names = (
    point_focal_source_bands
    + point_focal_band_names
    + point_alpha_output_bands
    + point_terrain_band_names
)
if len(point_expected_band_names) != 119:
    raise ValueError(
        f"Expected 119 predictors, found {len(point_expected_band_names)}"
    )
print("Canada-wide lazy sampling stack ready: 119 predictors + 5 QA bands")


# %% --------------------------------------------------------------------------
# 12. SAMPLE AND EXPORT SEPARATE CSV FILES BY ECOSYSTEM GROUP
# -----------------------------------------------------------------------------

GCS_BUCKET = "gee_exports_cd"
GCS_PREFIX = f"{CONFIG['export_folder']}/CanadaCoreCovariatesByEcosystem"

# This controls table exports only. The earlier external-core tasks have
# already been started, so leave False unless they must be submitted again.
RUN_POINT_EXPORTS = False

# Bounded batches avoid sending one very large client-side FeatureCollection.
BATCH_SIZE = 500


def _clean_property(value):
    """Convert pandas/numpy scalars to JSON-safe Python values."""
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def _frame_to_feature_collection(frame):
    features = []
    for _, row in frame.iterrows():
        properties = {}
        for key, value in row.to_dict().items():
            clean_value = _clean_property(value)
            if clean_value is not None:
                properties[key] = clean_value
        geometry = ee.Geometry.Point(
            [float(row["longitude"]), float(row["latitude"])]
        )
        features.append(ee.Feature(geometry, properties))
    return ee.FeatureCollection(features)


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")


core_export_tasks = []
if RUN_POINT_EXPORTS:
    for ecosystem_group, ecosystem_frame in selected_cores_df.groupby(
        "broad_ecosystem_group", dropna=False
    ):
        ecosystem_name = _safe_name(ecosystem_group or "unclassified")
        ecosystem_frame = ecosystem_frame.reset_index(drop=True)
        for batch_start in range(0, len(ecosystem_frame), BATCH_SIZE):
            batch_number = batch_start // BATCH_SIZE + 1
            batch = ecosystem_frame.iloc[
                batch_start : batch_start + BATCH_SIZE
            ].copy()
            batch_points = _frame_to_feature_collection(batch)
            batch_samples = point_sampling_image.sampleRegions(
                collection=batch_points,
                projection=FOCAL_PROJECTION,
                tileScale=8,
                geometries=True,
            )
            prefix = (
                f"RouteL_CanadaCoreCovariates_{ecosystem_name}_"
                f"batch{batch_number:02d}"
            )
            task = ee.batch.Export.table.toCloudStorage(
                collection=batch_samples,
                description=f"{prefix}_toGCS",
                bucket=GCS_BUCKET,
                fileNamePrefix=f"{GCS_PREFIX}/{prefix}",
                fileFormat="CSV",
            )
            task.start()
            core_export_tasks.append(
                {
                    "ecosystem_group": ecosystem_group,
                    "batch": batch_number,
                    "expected_rows": len(batch),
                    "prefix": prefix,
                    "task": task,
                }
            )
            print(
                f"Started {task.id}: {ecosystem_name} batch {batch_number:02d} "
                f"({len(batch)} cores)"
            )

    print(
        f"Submitted {len(core_export_tasks)} ecosystem-specific table tasks"
    )
    print(
        "The completed files together form the selected Canadian-core table."
    )
else:
    print("Point export skipped: set RUN_POINT_EXPORTS = True when ready.")


# %% --------------------------------------------------------------------------
# 13. CHECK TASK STATES WITHOUT RESTARTING EXPORTS
# -----------------------------------------------------------------------------

if not core_export_tasks:
    print("No point-export tasks were submitted in this run.")
else:
    for item in core_export_tasks:
        task_status = item["task"].status()
        print(
            item["ecosystem_group"],
            f"batch {item['batch']:02d}",
            f"expected rows {item['expected_rows']}",
            task_status.get("state"),
            task_status.get("error_message"),
        )
