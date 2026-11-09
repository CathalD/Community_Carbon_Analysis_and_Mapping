"""Sample the hierarchical-similarity pilot stack at Canadian donor cores.

Run this file in Google Colab. It builds the same 47 default continuous/QA
bands named by ``hierarchical_similarity_micro_pilot.js``, but as a lazy
Canada-wide point-sampling image. Carbon-stock response, coverage, and
completeness columns are never included in the Earth Engine inputs or outputs.

The default first test selects Riley and CanPEAT, groups by database source and
broad ecosystem, and splits each group into bounded client-side batches. The
script always prints and writes an export plan before it can submit tasks or
make synchronous downloads. Leave ``RUN_POINT_EXPORTS = False`` until the plan
and server band validation have been reviewed. The default DIRECT destination
writes donor CSV shards plus eligible community validation points to Colab's
temporary ``/content`` disk, verifies them, zips them, and offers one browser
download. Drive and Cloud Storage batch modes remain available.

An optional direct AOI raster path makes Route-L-aligned 4.5 km GeoTIFF tiles,
because Earth Engine direct image requests are limited to small chunks. The
tiles and their manifest are added to the same downloadable ZIP.

Important parity note
---------------------
Both this sampler and the paired pilot JavaScript explicitly assign the Route L
30 m projection to Sentinel-2 and Dynamic World temporal composites before
their focal reductions. Earth Engine collection reductions can otherwise
expose a generic WGS84 default projection. Band names, dates, masks, reducers,
radii, and source collections are kept identical. Numeric parity still needs
to be checked on a few shared validation points after the first Earth Engine
run.

Deliberate exclusions
---------------------
* Raw TWI: earlier AOI testing found it unstable in low-relief terrain.
* Optional coast, river, parent-material and permafrost hooks: they are not in
  the JavaScript's 47-band default stack and require separate source QA.
* Carbon response fields: properties are transferred with a strict metadata
  allow-list, not by copying every input CSV column.
"""


# %% --------------------------------------------------------------------------
# 0. IMPORTS, AUTHENTICATION AND USER SETTINGS
# -----------------------------------------------------------------------------

import hashlib
import json
import math
import os
from pathlib import Path
import re
import time
import urllib.request
import zipfile

import ee
import pandas as pd


# Change this only if a different Earth Engine billing/quota project is used.
GEE_PROJECT = "northstarlabs"

try:
    ee.Initialize(project=GEE_PROJECT)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=GEE_PROJECT)

try:
    from google.colab import files

    IN_COLAB = True
except ImportError:
    files = None
    IN_COLAB = False


CONFIG = {
    "start_year": 2019,
    "end_year": 2023,
    "s2_cloud_percent": 40,
    "scale_m": 30,
    "local_context_radius_m": 120,
    "regional_context_radius_m": 1000,
    "water_distance_max_m": 12000,
    "missing_predictor_value": -9999,
    "drive_folder": "MultiScaleCarbonCovariates",
    "gcs_bucket": "gee_exports_cd",
    "gcs_prefix": (
        "MultiScaleCarbonCovariates/HierarchicalSimilarityDonorCovariates"
    ),
    "export_prefix": "HierarchicalSimilarity_DonorCovariates",
    # Direct downloads are written to Colab's temporary local disk. Download
    # the resulting ZIP before the Colab runtime is disconnected.
    "direct_output_dir": "/content/HierarchicalSimilarity_DirectDownloads",
    "direct_zip_name": "HierarchicalSimilarity_DirectDownloads.zip",
    # Optional direct raster tiles use the configured Fort Severn AOI. A
    # 4.5 km tile is 150 x 150 cells at 30 m: about 4.2 MB uncompressed for
    # 47 float32 bands, comfortably below getDownloadURL's 32 MB limit.
    "direct_raster_aoi_asset": (
        "projects/north-star-project-470316/assets/aoi_fs"
    ),
    "direct_raster_tile_size_m": 4500,
    "direct_raster_prefix": "HierarchicalSimilarity_FS_AOI_Raster_30m",
}

if CONFIG["end_year"] < CONFIG["start_year"]:
    raise ValueError("CONFIG['end_year'] must be >= CONFIG['start_year']")

# These defaults are the cross-language parity contract. Changing the radii or
# scale is a new stack version and requires regenerating the explicit manifest.
if (
    CONFIG["scale_m"] != 30
    or CONFIG["local_context_radius_m"] != 120
    or CONFIG["regional_context_radius_m"] != 1000
):
    raise ValueError(
        "The v1 47-band parity contract requires 30 m, 120 m and 1000 m."
    )


# Exact Route L 30 m grid (custom Canada Lambert Conformal Conic / ESRI:102002).
ROUTE_L_CRS = (
    'PROJCS["Canada_Lambert_Conformal_Conic",'
    'GEOGCS["GCS_North_American_1983",'
    'DATUM["D_North_American_1983",'
    'SPHEROID["GRS_1980",6378137,298.257222101]],'
    'PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]],'
    'PROJECTION["Lambert_Conformal_Conic_2SP"],'
    'PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],'
    'PARAMETER["Central_Meridian",-96],'
    'PARAMETER["Standard_Parallel_1",50],'
    'PARAMETER["Standard_Parallel_2",70],'
    'PARAMETER["Latitude_Of_Origin",40],UNIT["Meter",1]]'
)
ROUTE_L_TRANSFORM_30M = [
    30.0,
    0.0,
    462400.096580488,
    0.0,
    -30.0,
    1881190.18705772,
]
OUTPUT_PROJECTION = ee.Projection(ROUTE_L_CRS, ROUTE_L_TRANSFORM_30M)


# %% --------------------------------------------------------------------------
# 1. LOAD AND SAFELY SELECT DONOR METADATA
# -----------------------------------------------------------------------------

PROJECT_INPUT_PATH = (
    "data/derived/external_carbon_profiles/"
    "external_carbon_cores_canada_for_gee.csv"
)
COLAB_INPUT_PATH = "/content/external_carbon_cores_canada_for_gee.csv"

if os.path.exists(COLAB_INPUT_PATH):
    core_table_path = COLAB_INPUT_PATH
elif os.path.exists(PROJECT_INPUT_PATH):
    core_table_path = PROJECT_INPUT_PATH
elif IN_COLAB:
    print("Upload external_carbon_cores_canada_for_gee.csv")
    uploaded = files.upload()
    if "external_carbon_cores_canada_for_gee.csv" in uploaded:
        core_table_path = "external_carbon_cores_canada_for_gee.csv"
    elif len(uploaded) == 1:
        core_table_path = next(iter(uploaded))
    else:
        raise FileNotFoundError("Could not identify the donor-core CSV upload")
else:
    raise FileNotFoundError(
        "Donor-core CSV not found. Run from the project root, upload it in "
        "Colab, or change PROJECT_INPUT_PATH."
    )

required_input_columns = {
    "core_uid",
    "database_code",
    "study_id",
    "latitude",
    "longitude",
    "country",
    "source_ecosystem",
    "broad_ecosystem_group",
}

# This is the complete set of source-table properties allowed into an export.
# Latitude/longitude are retained for joins and checking; response variables are
# never copied from the input DataFrame into Earth Engine features.
METADATA_COLUMNS = [
    "core_uid",
    "database_code",
    "study_id",
    "source_core_id",
    "source_ecosystem",
    "broad_ecosystem_group",
    "country",
    "region",
    "sampling_year",
    "latitude",
    "longitude",
]

forbidden_property_pattern = re.compile(
    r"(^carbon_|stock|coverage|complete_interval|profile_|n_layers|quality)",
    flags=re.IGNORECASE,
)
unsafe_metadata_columns = [
    name for name in METADATA_COLUMNS if forbidden_property_pattern.search(name)
]
if unsafe_metadata_columns:
    raise RuntimeError(
        "Unsafe properties in METADATA_COLUMNS: "
        + ", ".join(unsafe_metadata_columns)
    )

cores_df = pd.read_csv(core_table_path, low_memory=False)
missing_input_columns = required_input_columns - set(cores_df.columns)
if missing_input_columns:
    raise ValueError(
        "Rebuild the donor table; required columns are missing: "
        f"{sorted(missing_input_columns)}"
    )

missing_metadata_columns = set(METADATA_COLUMNS) - set(cores_df.columns)
if missing_metadata_columns:
    raise ValueError(
        "Metadata allow-list columns are missing: "
        f"{sorted(missing_metadata_columns)}"
    )

cores_df["database_code"] = cores_df["database_code"].astype(str).str.upper()
cores_df = cores_df[
    cores_df["country"].astype(str).str.casefold().eq("canada")
].copy()
cores_df = cores_df.dropna(subset=["latitude", "longitude"])
cores_df = cores_df[
    cores_df["latitude"].between(-90, 90)
    & cores_df["longitude"].between(-180, 180)
].drop_duplicates("core_uid")

# Tractable first test. Set to None only after Riley+CanPEAT outputs pass QA.
DATABASE_FILTER = ["RILEY", "CANPEAT"]

# Leave as None for every ecosystem represented by the selected sources.
ECOSYSTEM_GROUP_FILTER = None

selected_cores_df = cores_df.copy()
if DATABASE_FILTER is not None:
    selected_cores_df = selected_cores_df[
        selected_cores_df["database_code"].isin(DATABASE_FILTER)
    ]
if ECOSYSTEM_GROUP_FILTER is not None:
    selected_cores_df = selected_cores_df[
        selected_cores_df["broad_ecosystem_group"].isin(
            ECOSYSTEM_GROUP_FILTER
        )
    ]
selected_cores_df = selected_cores_df.reset_index(drop=True)
selected_cores_df["record_role"] = "donor"

if selected_cores_df.empty:
    raise ValueError("The source/ecosystem filters selected no donor cores")

print(f"Loaded {len(cores_df):,} unique Canadian donor cores")
print(f"Selected {len(selected_cores_df):,} donor cores for this run")
print(selected_cores_df["database_code"].value_counts().sort_index())
print(
    "Exported source properties are restricted to: "
    + ", ".join(METADATA_COLUMNS)
)
print("Carbon response columns transferred: 0")


# %% --------------------------------------------------------------------------
# 2. DEFINE THE EXACT 47-BAND PARITY MANIFEST
# -----------------------------------------------------------------------------

EXPECTED_BAND_NAMES = [
    "local_s2_ndmi_spring_median",
    "local_s2_ndmi_spring_temporal_sd",
    "local_s2_ndvi_summer_median",
    "local_s2_ndvi_summer_temporal_sd",
    "local_s2_evi_summer_median",
    "local_s2_evi_summer_temporal_sd",
    "local_s2_ndvi_120m_mean",
    "local_s2_ndvi_120m_spatial_sd",
    "local_dw_tree_prob_mean",
    "local_dw_tree_prob_temporal_sd",
    "local_dw_tree_120m_mean",
    "local_dw_shrub_prob_mean",
    "local_dw_shrub_prob_temporal_sd",
    "local_dw_shrub_120m_mean",
    "local_dw_floodedveg_prob_mean",
    "local_dw_floodedveg_prob_temporal_sd",
    "local_dw_floodedveg_120m_mean",
    "local_dw_water_prob_mean",
    "local_dw_water_prob_temporal_sd",
    "local_dw_water_120m_mean",
    "local_s1_vv_db_median",
    "local_s1_vv_db_temporal_sd",
    "local_s1_vh_db_median",
    "local_s1_vh_db_temporal_sd",
    "local_s1_vv_minus_vh_db_median",
    "local_s1_vv_minus_vh_db_temporal_sd",
    "local_s1_vv_over_vh_linear_median",
    "local_s1_vv_over_vh_linear_temporal_sd",
    "local_slope_deg",
    "local_relative_elevation_120m",
    "local_elevation_120m_sd",
    "regional_jrc_water_occurrence_pct",
    "regional_jrc_water_seasonality_months",
    "regional_jrc_water_recurrence_pct",
    "regional_observed_water_distance_m",
    "regional_elevation_m",
    "regional_relative_elevation_1000m",
    "regional_elevation_1000m_sd",
    "broad_era5_mean_temperature_c",
    "broad_era5_monthly_temperature_sd_c",
    "broad_era5_mean_annual_precip_mm",
    "broad_era5_annual_precip_sd_mm",
    "broad_era5_mean_freezing_months",
    "qa_s2_spring_observations",
    "qa_s2_summer_observations",
    "qa_dynamic_world_observations",
    "qa_sentinel1_observations",
]

if len(EXPECTED_BAND_NAMES) != 47:
    raise RuntimeError(
        f"Expected a 47-band manifest, found {len(EXPECTED_BAND_NAMES)}"
    )
if len(EXPECTED_BAND_NAMES) != len(set(EXPECTED_BAND_NAMES)):
    raise RuntimeError("Duplicate names found in the 47-band manifest")
if not all(
    name.startswith(("local_", "regional_", "broad_", "qa_"))
    for name in EXPECTED_BAND_NAMES
):
    raise RuntimeError("Every band must have a scale-family or QA prefix")

band_manifest_df = pd.DataFrame(
    {
        "band_order": range(1, len(EXPECTED_BAND_NAMES) + 1),
        "band_name": EXPECTED_BAND_NAMES,
        "block": [name.split("_", 1)[0] for name in EXPECTED_BAND_NAMES],
        "response_field": False,
        "definition_version": "hierarchical_similarity_micro_pilot_v1",
        "parity_note": (
            "same intended definition; both exporters pin 30m projection "
            "before focal reduction"
        ),
    }
)

EXPORT_METADATA_COLUMNS = METADATA_COLUMNS + ["record_role"]
planned_export_fields = EXPORT_METADATA_COLUMNS + EXPECTED_BAND_NAMES
unsafe_planned_fields = [
    name
    for name in planned_export_fields
    if forbidden_property_pattern.search(name)
]
if unsafe_planned_fields:
    raise RuntimeError(
        "A response-like field would enter the export: "
        + ", ".join(unsafe_planned_fields)
    )
BAND_MANIFEST_PATH = "/content/hierarchical_similarity_band_manifest_v1.csv"
if not IN_COLAB:
    BAND_MANIFEST_PATH = "/tmp/hierarchical_similarity_band_manifest_v1.csv"
band_manifest_df.to_csv(BAND_MANIFEST_PATH, index=False)
print(f"Wrote 47-band manifest: {BAND_MANIFEST_PATH}")


# %% --------------------------------------------------------------------------
# 3. EARTH ENGINE HELPERS
# -----------------------------------------------------------------------------

CANADA_BOUNDS = ee.Geometry.Rectangle(
    [-141.1, 41.0, -52.5, 83.2], proj="EPSG:4326", geodesic=False
)
START_DATE = ee.Date.fromYMD(CONFIG["start_year"], 1, 1)
END_DATE = ee.Date.fromYMD(CONFIG["end_year"] + 1, 1, 1)


def _combine_images(images):
    """Concatenate a non-empty Python list of images without name prefixes."""
    if not images:
        raise ValueError("At least one image is required")
    result = ee.Image(images[0])
    for image in images[1:]:
        result = result.addBands(ee.Image(image))
    return result


def _mask_sentinel2(image):
    """Keep SCL vegetation, bare-ground and water pixels; scale reflectance."""
    scl = image.select("SCL")
    clear = scl.eq(4).Or(scl.eq(5)).Or(scl.eq(6))
    return (
        image.select(["B2", "B3", "B4", "B8", "B11"])
        .multiply(0.0001)
        .updateMask(clear)
        .copyProperties(image, ["system:time_start"])
    )


def _add_sentinel2_indices(image):
    ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI")
    ndmi = image.normalizedDifference(["B8", "B11"]).rename("NDMI")
    evi = image.expression(
        "2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)",
        {
            "nir": image.select("B8"),
            "red": image.select("B4"),
            "blue": image.select("B2"),
        },
    ).rename("EVI")
    return image.addBands([ndvi, ndmi, evi])


def _focal_mean(image, radius_m, output_name):
    return image.reduceNeighborhood(
        reducer=ee.Reducer.mean(),
        kernel=ee.Kernel.circle(radius=radius_m, units="meters"),
    ).rename(output_name)


def _focal_std_dev(image, radius_m, output_name):
    return image.reduceNeighborhood(
        reducer=ee.Reducer.stdDev(),
        kernel=ee.Kernel.circle(radius=radius_m, units="meters"),
    ).rename(output_name)


def _add_sentinel1_metrics(image):
    contrast_db = image.select("VV").subtract(image.select("VH")).rename(
        "VV_minus_VH_dB"
    )
    ratio_linear = ee.Image.constant(10).pow(contrast_db.divide(10)).rename(
        "VV_over_VH_linear"
    )
    return (
        image.select(["VV", "VH"])
        .addBands([contrast_db, ratio_linear])
        .copyProperties(image, ["system:time_start"])
    )


# %% --------------------------------------------------------------------------
# 4. LOCAL BLOCK: S2, DYNAMIC WORLD, S1 AND TERRAIN
# -----------------------------------------------------------------------------

s2 = (
    ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(CANADA_BOUNDS)
    .filterDate(START_DATE, END_DATE)
    .filter(
        ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", CONFIG["s2_cloud_percent"])
    )
    .map(_mask_sentinel2)
    .map(_add_sentinel2_indices)
)
s2_spring = s2.filter(ee.Filter.calendarRange(4, 6, "month"))
s2_summer = s2.filter(ee.Filter.calendarRange(6, 8, "month"))

# Explicit projection pinning is the operational difference documented above.
ndmi_spring_median = (
    s2_spring.select("NDMI").median().setDefaultProjection(OUTPUT_PROJECTION)
)
ndvi_summer_median = (
    s2_summer.select("NDVI").median().setDefaultProjection(OUTPUT_PROJECTION)
)
evi_summer_median = (
    s2_summer.select("EVI").median().setDefaultProjection(OUTPUT_PROJECTION)
)

local_images = [
    ndmi_spring_median.rename("local_s2_ndmi_spring_median"),
    s2_spring.select("NDMI")
    .reduce(ee.Reducer.stdDev())
    .rename("local_s2_ndmi_spring_temporal_sd"),
    ndvi_summer_median.rename("local_s2_ndvi_summer_median"),
    s2_summer.select("NDVI")
    .reduce(ee.Reducer.stdDev())
    .rename("local_s2_ndvi_summer_temporal_sd"),
    evi_summer_median.rename("local_s2_evi_summer_median"),
    s2_summer.select("EVI")
    .reduce(ee.Reducer.stdDev())
    .rename("local_s2_evi_summer_temporal_sd"),
    _focal_mean(
        ndvi_summer_median,
        CONFIG["local_context_radius_m"],
        "local_s2_ndvi_120m_mean",
    ),
    _focal_std_dev(
        ndvi_summer_median,
        CONFIG["local_context_radius_m"],
        "local_s2_ndvi_120m_spatial_sd",
    ),
]

dynamic_world = (
    ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
    .filterBounds(CANADA_BOUNDS)
    .filterDate(START_DATE, END_DATE)
    .filter(ee.Filter.calendarRange(6, 8, "month"))
)
dw_bands = ["trees", "shrub_and_scrub", "flooded_vegetation", "water"]
dw_names = ["tree", "shrub", "floodedveg", "water"]
for source_band, output_root in zip(dw_bands, dw_names):
    dw_mean = (
        dynamic_world.select(source_band)
        .mean()
        .setDefaultProjection(OUTPUT_PROJECTION)
    )
    local_images.extend(
        [
            dw_mean.rename(f"local_dw_{output_root}_prob_mean"),
            dynamic_world.select(source_band)
            .reduce(ee.Reducer.stdDev())
            .rename(f"local_dw_{output_root}_prob_temporal_sd"),
            _focal_mean(
                dw_mean,
                CONFIG["local_context_radius_m"],
                f"local_dw_{output_root}_120m_mean",
            ),
        ]
    )

sentinel1 = (
    ee.ImageCollection("COPERNICUS/S1_GRD")
    .filterBounds(CANADA_BOUNDS)
    .filterDate(START_DATE, END_DATE)
    .filter(ee.Filter.eq("instrumentMode", "IW"))
    .filter(
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV")
    )
    .filter(
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH")
    )
    .map(_add_sentinel1_metrics)
)
s1_metric_bands = ["VV", "VH", "VV_minus_VH_dB", "VV_over_VH_linear"]
s1_metric_names = [
    "vv_db",
    "vh_db",
    "vv_minus_vh_db",
    "vv_over_vh_linear",
]
for source_band, output_root in zip(s1_metric_bands, s1_metric_names):
    local_images.extend(
        [
            sentinel1.select(source_band)
            .median()
            .rename(f"local_s1_{output_root}_median"),
            sentinel1.select(source_band)
            .reduce(ee.Reducer.stdDev())
            .rename(f"local_s1_{output_root}_temporal_sd"),
        ]
    )

cdem_collection = (
    ee.ImageCollection("NRCan/CDEM")
    .filterBounds(CANADA_BOUNDS)
    .select("elevation")
)
cdem_first = ee.Image(cdem_collection.first())
elevation = (
    cdem_collection.mosaic()
    .setDefaultProjection(cdem_first.projection())
    .resample("bilinear")
    .rename("elevation_m")
)
slope = ee.Terrain.slope(elevation)
local_elevation_mean = _focal_mean(
    elevation,
    CONFIG["local_context_radius_m"],
    "local_elevation_context_mean",
)
local_images.extend(
    [
        slope.rename("local_slope_deg"),
        elevation.subtract(local_elevation_mean).rename(
            "local_relative_elevation_120m"
        ),
        _focal_std_dev(
            elevation,
            CONFIG["local_context_radius_m"],
            "local_elevation_120m_sd",
        ),
    ]
)
local_stack = _combine_images(local_images).toFloat()


# %% --------------------------------------------------------------------------
# 5. REGIONAL BLOCK: WATER REGIME AND MESOTOPOGRAPHY
# -----------------------------------------------------------------------------

surface_water = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
# JRC masks locations where surface water was never detected. For these three
# occurrence/frequency descriptors that state is a meaningful zero rather than
# an unknown value. Keep this identical to the micro-pilot raster definition.
regional_images = [
    surface_water.select("occurrence").unmask(0).rename(
        "regional_jrc_water_occurrence_pct"
    ),
    surface_water.select("seasonality").unmask(0).rename(
        "regional_jrc_water_seasonality_months"
    ),
    surface_water.select("recurrence").unmask(0).rename(
        "regional_jrc_water_recurrence_pct"
    ),
]

max_water_distance_pixels = math.ceil(
    CONFIG["water_distance_max_m"] / CONFIG["scale_m"]
)
water_binary_30m = (
    surface_water.select("max_extent")
    .eq(1)
    .unmask(0)
    .toByte()
    .setDefaultProjection(OUTPUT_PROJECTION)
)
water_distance = (
    water_binary_30m.fastDistanceTransform(
        max_water_distance_pixels,
        "pixels",
        "squared_euclidean",
    )
    .sqrt()
    .multiply(CONFIG["scale_m"])
    .min(CONFIG["water_distance_max_m"])
    .rename("regional_observed_water_distance_m")
)
regional_images.append(water_distance)

regional_elevation_mean = _focal_mean(
    elevation,
    CONFIG["regional_context_radius_m"],
    "regional_elevation_context_mean",
)
regional_images.extend(
    [
        elevation.rename("regional_elevation_m"),
        elevation.subtract(regional_elevation_mean).rename(
            "regional_relative_elevation_1000m"
        ),
        _focal_std_dev(
            elevation,
            CONFIG["regional_context_radius_m"],
            "regional_elevation_1000m_sd",
        ),
    ]
)
regional_stack = _combine_images(regional_images).toFloat()


# %% --------------------------------------------------------------------------
# 6. BROAD BLOCK: MULTI-YEAR ERA5-LAND CLIMATE
# -----------------------------------------------------------------------------

era5 = (
    ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")
    .filterBounds(CANADA_BOUNDS)
    .filterDate(START_DATE, END_DATE)
)
years = ee.List.sequence(CONFIG["start_year"], CONFIG["end_year"])


def _annual_precipitation(year):
    year = ee.Number(year)
    year_start = ee.Date.fromYMD(year, 1, 1)
    return (
        era5.filterDate(year_start, year_start.advance(1, "year"))
        .select("total_precipitation_sum")
        .sum()
        .multiply(1000)
        .rename("annual_precip_mm")
        .set("year", year)
    )


def _annual_freezing_months(year):
    year = ee.Number(year)
    year_start = ee.Date.fromYMD(year, 1, 1)
    return (
        era5.filterDate(year_start, year_start.advance(1, "year"))
        .select("temperature_2m")
        .map(lambda image: image.lt(273.15))
        .sum()
        .rename("freezing_months")
        .set("year", year)
    )


annual_precipitation = ee.ImageCollection.fromImages(
    years.map(_annual_precipitation)
)
annual_freezing_months = ee.ImageCollection.fromImages(
    years.map(_annual_freezing_months)
)

broad_stack = _combine_images(
    [
        era5.select("temperature_2m")
        .mean()
        .subtract(273.15)
        .rename("broad_era5_mean_temperature_c"),
        era5.select("temperature_2m")
        .reduce(ee.Reducer.stdDev())
        .rename("broad_era5_monthly_temperature_sd_c"),
        annual_precipitation.mean().rename(
            "broad_era5_mean_annual_precip_mm"
        ),
        annual_precipitation.reduce(ee.Reducer.stdDev()).rename(
            "broad_era5_annual_precip_sd_mm"
        ),
        annual_freezing_months.mean().rename(
            "broad_era5_mean_freezing_months"
        ),
    ]
).toFloat()


# %% --------------------------------------------------------------------------
# 7. QA, FINAL ORDER AND SERVER-SIDE BAND VALIDATION
# -----------------------------------------------------------------------------

qa_stack = (
    _combine_images(
        [
            s2_spring.select("NDMI")
            .count()
            .rename("qa_s2_spring_observations"),
            s2_summer.select("NDVI")
            .count()
            .rename("qa_s2_summer_observations"),
            dynamic_world.select("trees")
            .count()
            .rename("qa_dynamic_world_observations"),
            sentinel1.select("VV")
            .count()
            .rename("qa_sentinel1_observations"),
        ]
    )
    .unmask(0)
    .toFloat()
)

covariate_stack = (
    local_stack.addBands(regional_stack)
    .addBands(broad_stack)
    .addBands(qa_stack)
    .select(EXPECTED_BAND_NAMES)
    .set(
        {
            "purpose": "hierarchical multi-scale donor point sampling",
            "definition_version": "hierarchical_similarity_micro_pilot_v1",
            "scale_families": "local|regional|broad|qa",
            "start_year": CONFIG["start_year"],
            "end_year": CONFIG["end_year"],
            "output_scale_m": CONFIG["scale_m"],
            "raw_twi_included": False,
            "carbon_response_fields_included": False,
            "projection_parity_note": (
                "30m default projection pinned before focal reductions"
            ),
        }
    )
)

# This read-only metadata request does not sample pixels or submit exports.
VALIDATE_SERVER_BANDS = True
if VALIDATE_SERVER_BANDS:
    server_band_names = covariate_stack.bandNames().getInfo()
    if server_band_names != EXPECTED_BAND_NAMES:
        missing = sorted(set(EXPECTED_BAND_NAMES) - set(server_band_names))
        unexpected = sorted(set(server_band_names) - set(EXPECTED_BAND_NAMES))
        raise RuntimeError(
            "Server stack does not match the 47-band manifest. "
            f"Missing={missing}; unexpected={unexpected}"
        )
    print("Server band validation passed: exact 47-band names and order")
else:
    print("WARNING: server band validation disabled")

# Preserve every donor row even if a predictor is masked. QA counts distinguish
# weak temporal support; -9999 must be converted to missing during modelling.
qa_band_names = [
    name for name in EXPECTED_BAND_NAMES if name.startswith("qa_")
]
predictor_band_names = [
    name for name in EXPECTED_BAND_NAMES if not name.startswith("qa_")
]
point_sampling_image = covariate_stack.select(predictor_band_names).unmask(
    CONFIG["missing_predictor_value"]
).addBands(covariate_stack.select(qa_band_names).unmask(0))


# %% --------------------------------------------------------------------------
# 8. BUILD THE SOURCE x ECOSYSTEM x BATCH EXPORT PLAN
# -----------------------------------------------------------------------------

BATCH_SIZE = 250
if BATCH_SIZE < 1:
    raise ValueError("BATCH_SIZE must be positive")


def _safe_name(value):
    if pd.isna(value) or str(value).strip() == "":
        return "unclassified"
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")


plan_rows = []
group_columns = ["database_code", "broad_ecosystem_group"]
for group_keys, group_frame in selected_cores_df.groupby(
    group_columns, dropna=False, sort=True
):
    database_code, ecosystem_group = group_keys
    n_cores = len(group_frame)
    n_batches = math.ceil(n_cores / BATCH_SIZE)
    for batch_number in range(1, n_batches + 1):
        row_start = (batch_number - 1) * BATCH_SIZE
        expected_rows = min(BATCH_SIZE, n_cores - row_start)
        plan_rows.append(
            {
                "database_code": database_code,
                "broad_ecosystem_group": ecosystem_group,
                "batch": batch_number,
                "expected_rows": expected_rows,
                "n_group_cores": n_cores,
                "n_group_batches": n_batches,
            }
        )

export_plan_df = pd.DataFrame(plan_rows)
EXPORT_PLAN_PATH = "/content/hierarchical_similarity_donor_export_plan.csv"
if not IN_COLAB:
    EXPORT_PLAN_PATH = "/tmp/hierarchical_similarity_donor_export_plan.csv"
export_plan_df.to_csv(EXPORT_PLAN_PATH, index=False)

print("\nEXPORT PLAN — no tasks have been submitted")
print(export_plan_df.to_string(index=False))
print(f"Planned files: {len(export_plan_df):,}")
print(f"Planned donor rows: {export_plan_df['expected_rows'].sum():,}")
print(f"Wrote export plan: {EXPORT_PLAN_PATH}")

if export_plan_df["expected_rows"].sum() != len(selected_cores_df):
    raise RuntimeError("Export-plan rows do not sum to the selected donor count")


# %% --------------------------------------------------------------------------
# 9. SAMPLE AND EXPORT — DISABLED BY DEFAULT
# -----------------------------------------------------------------------------

# DIRECT writes small CSV shards to /content and makes one browser-downloadable
# ZIP. DRIVE and GCS submit normal Earth Engine batch tasks. Direct mode avoids
# both services, but Colab's /content disk is temporary.
EXPORT_DESTINATION = "DIRECT"
if EXPORT_DESTINATION not in {"DIRECT", "DRIVE", "GCS"}:
    raise ValueError("EXPORT_DESTINATION must be 'DIRECT', 'DRIVE' or 'GCS'")

# Review Sections 1, 2, 7 and 8 before enabling either switch. The raster
# option is separate because a 47-band AOI raster must be downloaded in tiles.
RUN_POINT_EXPORTS = False
INCLUDE_COMMUNITY_VALIDATION_POINTS = True
RUN_DIRECT_RASTER_TILES = False
DOWNLOAD_DIRECT_ZIP_TO_BROWSER = True

if RUN_DIRECT_RASTER_TILES and EXPORT_DESTINATION != "DIRECT":
    raise ValueError("Raster tile downloads require EXPORT_DESTINATION='DIRECT'")


def _clean_property(value):
    """Convert pandas/numpy scalars to Earth Engine-safe values."""
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def _frame_to_feature_collection(frame):
    """Create points with allow-listed metadata and no response properties."""
    features = []
    for _, row in frame.iterrows():
        properties = {}
        for key in EXPORT_METADATA_COLUMNS:
            clean_value = _clean_property(row[key])
            if clean_value is not None:
                properties[key] = clean_value
        point = ee.Geometry.Point(
            [float(row["longitude"]), float(row["latitude"])]
        )
        features.append(ee.Feature(point, properties))
    return ee.FeatureCollection(features)


def _download_url(url, destination, attempts=4):
    """Download one signed Earth Engine URL with bounded retries."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "hierarchical-similarity-colab/1.0"},
            )
            with urllib.request.urlopen(request, timeout=600) as response:
                with partial.open("wb") as output:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)
            partial.replace(destination)
            return destination
        except Exception as error:
            last_error = error
            if partial.exists():
                partial.unlink()
            if attempt < attempts:
                wait_seconds = 2 ** attempt
                print(
                    f"Download attempt {attempt} failed; retrying in "
                    f"{wait_seconds}s: {error}"
                )
                time.sleep(wait_seconds)
    raise RuntimeError(
        f"Download failed after {attempts} attempts: {destination}"
    ) from last_error


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _flatten_coordinate_pairs(value):
    """Return every numeric [x, y] pair in a nested EE coordinate list."""
    if (
        isinstance(value, list)
        and len(value) >= 2
        and all(isinstance(item, (int, float)) for item in value[:2])
    ):
        return [(float(value[0]), float(value[1]))]
    pairs = []
    if isinstance(value, list):
        for item in value:
            pairs.extend(_flatten_coordinate_pairs(item))
    return pairs


def _is_geotiff(path):
    with Path(path).open("rb") as handle:
        signature = handle.read(4)
    return signature in {b"II*\x00", b"MM\x00*"}


direct_root = None
direct_point_manifest = []
direct_raster_manifest = []
if EXPORT_DESTINATION == "DIRECT" and (
    RUN_POINT_EXPORTS or RUN_DIRECT_RASTER_TILES
):
    run_id = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    configured_root = Path(CONFIG["direct_output_dir"])
    if not IN_COLAB:
        configured_root = Path("/tmp/HierarchicalSimilarity_DirectDownloads")
    direct_root = configured_root.with_name(
        configured_root.name + "_" + run_id
    )
    direct_root.mkdir(parents=True, exist_ok=False)
    band_manifest_df.to_csv(direct_root / "band_manifest.csv", index=False)
    export_plan_df.to_csv(direct_root / "point_export_plan.csv", index=False)
    print(f"Direct-download run directory: {direct_root}")


export_tasks = []
if RUN_POINT_EXPORTS:
    if EXPORT_DESTINATION == "DIRECT":
        print(f"Downloading {len(export_plan_df)} CSV shards to Colab...")
    else:
        print(
            f"Submitting {len(export_plan_df)} {EXPORT_DESTINATION} "
            "table tasks..."
        )
    for group_keys, group_frame in selected_cores_df.groupby(
        group_columns, dropna=False, sort=True
    ):
        database_code, ecosystem_group = group_keys
        database_name = _safe_name(database_code)
        ecosystem_name = _safe_name(ecosystem_group)
        group_frame = group_frame.reset_index(drop=True)

        for batch_start in range(0, len(group_frame), BATCH_SIZE):
            batch_number = batch_start // BATCH_SIZE + 1
            batch = group_frame.iloc[
                batch_start : batch_start + BATCH_SIZE
            ].copy()
            batch_points = _frame_to_feature_collection(batch)
            batch_samples = point_sampling_image.sampleRegions(
                collection=batch_points,
                projection=OUTPUT_PROJECTION,
                tileScale=8,
                geometries=False,
            )
            prefix = (
                f"{CONFIG['export_prefix']}_{database_name}_"
                f"{ecosystem_name}_batch{batch_number:02d}"
            )

            if EXPORT_DESTINATION == "DIRECT":
                destination = direct_root / "point_shards" / f"{prefix}.csv"
                download_url = batch_samples.getDownloadURL(
                    filetype="CSV",
                    selectors=planned_export_fields,
                    filename=prefix,
                )
                _download_url(download_url, destination)
                downloaded = pd.read_csv(destination, low_memory=False)
                missing_fields = sorted(
                    set(planned_export_fields) - set(downloaded.columns)
                )
                if missing_fields:
                    raise RuntimeError(
                        f"{destination.name} is missing fields: {missing_fields}"
                    )
                if len(downloaded) != len(batch):
                    raise RuntimeError(
                        f"{destination.name} has {len(downloaded)} rows; "
                        f"expected {len(batch)}"
                    )
                direct_point_manifest.append(
                    {
                        "database_code": str(database_code),
                        "broad_ecosystem_group": str(ecosystem_group),
                        "batch": int(batch_number),
                        "rows": int(len(downloaded)),
                        "file": str(destination.relative_to(direct_root)),
                        "bytes": int(destination.stat().st_size),
                        "sha256": _sha256(destination),
                    }
                )
                print(
                    f"Downloaded {destination.name}: "
                    f"{len(downloaded)} rows"
                )
                continue

            if EXPORT_DESTINATION == "DRIVE":
                task = ee.batch.Export.table.toDrive(
                    collection=batch_samples,
                    description=f"{prefix}_toDrive",
                    folder=CONFIG["drive_folder"],
                    fileNamePrefix=prefix,
                    fileFormat="CSV",
                    selectors=planned_export_fields,
                )
            else:
                task = ee.batch.Export.table.toCloudStorage(
                    collection=batch_samples,
                    description=f"{prefix}_toGCS",
                    bucket=CONFIG["gcs_bucket"],
                    fileNamePrefix=f"{CONFIG['gcs_prefix']}/{prefix}",
                    fileFormat="CSV",
                    selectors=planned_export_fields,
                )

            task.start()
            export_tasks.append(
                {
                    "database_code": database_code,
                    "broad_ecosystem_group": ecosystem_group,
                    "batch": batch_number,
                    "expected_rows": len(batch),
                    "prefix": prefix,
                    "task": task,
                }
            )
            print(
                f"Started {task.id}: {database_name} / {ecosystem_name} / "
                f"batch {batch_number:02d} ({len(batch)} rows)"
            )
    if EXPORT_DESTINATION == "DIRECT":
        print(f"Downloaded {len(direct_point_manifest)} point-table shards")
    else:
        print(f"Submitted {len(export_tasks)} table tasks")
else:
    print("Point exports skipped: RUN_POINT_EXPORTS remains False")


# Add the community points to the same direct-download package with an explicit
# validation role. Carbon measurements are deliberately not uploaded to Earth
# Engine; the local R step joins them later by core_uid.
if RUN_POINT_EXPORTS and INCLUDE_COMMUNITY_VALIDATION_POINTS:
    if EXPORT_DESTINATION != "DIRECT":
        raise ValueError(
            "Community validation points are packaged by the DIRECT workflow. "
            "Set INCLUDE_COMMUNITY_VALIDATION_POINTS=False for batch exports."
        )
    community_candidates = [
        Path("/content/community_cores_for_gee.csv"),
        Path("data/derived/community_cores_for_gee.csv"),
    ]
    community_path = next(
        (path for path in community_candidates if path.exists()), None
    )
    if community_path is None and IN_COLAB:
        print("Upload community_cores_for_gee.csv")
        uploaded = files.upload()
        if "community_cores_for_gee.csv" in uploaded:
            community_path = Path("community_cores_for_gee.csv")
        elif len(uploaded) == 1:
            community_path = Path(next(iter(uploaded)))
    if community_path is None:
        raise FileNotFoundError(
            "community_cores_for_gee.csv not found. Upload it or set "
            "INCLUDE_COMMUNITY_VALIDATION_POINTS=False."
        )

    community_df = pd.read_csv(community_path, low_memory=False)
    required_community_columns = {
        "core_uid",
        "database_code",
        "source_core_id",
        "plot_id",
        "latitude",
        "longitude",
        "country",
        "region",
        "sampling_period",
        "source_ecosystem",
        "broad_ecosystem_group",
        "gee_eligible",
    }
    missing_community_columns = required_community_columns - set(
        community_df.columns
    )
    if missing_community_columns:
        raise ValueError(
            "Community point table is missing: "
            f"{sorted(missing_community_columns)}"
        )
    eligible = community_df["gee_eligible"].astype(str).str.casefold().isin(
        {"true", "t", "yes", "y", "1"}
    )
    community_df = community_df[eligible].copy()
    community_df = community_df.dropna(subset=["latitude", "longitude"])
    community_df = community_df[
        community_df["latitude"].between(-90, 90)
        & community_df["longitude"].between(-180, 180)
    ].drop_duplicates("core_uid")
    if community_df.empty:
        raise ValueError("No eligible community validation points remain")

    community_df["study_id"] = community_df["plot_id"]
    community_df["sampling_year"] = community_df["sampling_period"]
    community_df["record_role"] = "validation"
    missing_export_metadata = set(EXPORT_METADATA_COLUMNS) - set(
        community_df.columns
    )
    if missing_export_metadata:
        raise RuntimeError(
            "Community metadata could not be standardized: "
            f"{sorted(missing_export_metadata)}"
        )

    community_points = _frame_to_feature_collection(community_df)
    community_samples = point_sampling_image.sampleRegions(
        collection=community_points,
        projection=OUTPUT_PROJECTION,
        tileScale=8,
        geometries=False,
    )
    community_prefix = (
        f"{CONFIG['export_prefix']}_COMMUNITY_validation_batch01"
    )
    community_destination = (
        direct_root / "point_shards" / f"{community_prefix}.csv"
    )
    community_url = community_samples.getDownloadURL(
        filetype="CSV",
        selectors=planned_export_fields,
        filename=community_prefix,
    )
    _download_url(community_url, community_destination)
    community_downloaded = pd.read_csv(
        community_destination, low_memory=False
    )
    if len(community_downloaded) != len(community_df):
        raise RuntimeError(
            f"Community shard has {len(community_downloaded)} rows; "
            f"expected {len(community_df)}"
        )
    missing_fields = sorted(
        set(planned_export_fields) - set(community_downloaded.columns)
    )
    if missing_fields:
        raise RuntimeError(
            f"Community shard is missing fields: {missing_fields}"
        )
    direct_point_manifest.append(
        {
            "database_code": "COMMUNITY",
            "broad_ecosystem_group": "community_validation",
            "batch": 1,
            "rows": int(len(community_downloaded)),
            "file": str(community_destination.relative_to(direct_root)),
            "bytes": int(community_destination.stat().st_size),
            "sha256": _sha256(community_destination),
        }
    )
    print(
        f"Downloaded community validation shard: "
        f"{len(community_downloaded)} rows"
    )


# %% --------------------------------------------------------------------------
# 10. OPTIONAL DIRECT AOI RASTER DOWNLOAD IN SAFE, ALIGNED TILES
# -----------------------------------------------------------------------------

if RUN_DIRECT_RASTER_TILES:
    tile_size_m = int(CONFIG["direct_raster_tile_size_m"])
    scale_m = int(CONFIG["scale_m"])
    if tile_size_m <= 0 or tile_size_m % scale_m != 0:
        raise ValueError(
            "direct_raster_tile_size_m must be positive and divisible by 30"
        )

    tile_pixels = tile_size_m // scale_m
    uncompressed_tile_mb = (
        tile_pixels * tile_pixels * len(EXPECTED_BAND_NAMES) * 4
    ) / (1024 ** 2)
    if uncompressed_tile_mb >= 30:
        raise ValueError(
            "Direct tile is too close to Earth Engine's 32 MB request limit: "
            f"{uncompressed_tile_mb:.1f} MB"
        )

    raster_aoi = ee.FeatureCollection(
        CONFIG["direct_raster_aoi_asset"]
    ).geometry(1)
    projected_bounds = raster_aoi.bounds(1, ROUTE_L_CRS)
    pairs = _flatten_coordinate_pairs(projected_bounds.coordinates().getInfo())
    if not pairs:
        raise RuntimeError("Could not resolve projected AOI bounds")
    x_values = [pair[0] for pair in pairs]
    y_values = [pair[1] for pair in pairs]
    min_x, max_x = min(x_values), max(x_values)
    min_y, max_y = min(y_values), max(y_values)

    origin_x = float(ROUTE_L_TRANSFORM_30M[2])
    origin_y = float(ROUTE_L_TRANSFORM_30M[5])
    first_col = math.floor((min_x - origin_x) / tile_size_m)
    last_col = math.ceil((max_x - origin_x) / tile_size_m)
    first_row = math.floor((origin_y - max_y) / tile_size_m)
    last_row = math.ceil((origin_y - min_y) / tile_size_m)

    candidate_tiles = (
        max(0, last_col - first_col) * max(0, last_row - first_row)
    )
    print(
        f"Direct raster grid: up to {candidate_tiles} aligned tiles; "
        f"{tile_pixels} x {tile_pixels} cells; "
        f"~{uncompressed_tile_mb:.1f} MB uncompressed per tile"
    )

    raster_image = point_sampling_image.clip(raster_aoi).toFloat()
    raster_directory = direct_root / "raster_tiles"
    raster_directory.mkdir(parents=True, exist_ok=True)

    for row_index in range(first_row, last_row):
        tile_ymax = origin_y - row_index * tile_size_m
        tile_ymin = tile_ymax - tile_size_m
        for col_index in range(first_col, last_col):
            tile_xmin = origin_x + col_index * tile_size_m
            tile_xmax = tile_xmin + tile_size_m
            tile_geometry = ee.Geometry.Rectangle(
                [tile_xmin, tile_ymin, tile_xmax, tile_ymax],
                proj=ROUTE_L_CRS,
                geodesic=False,
            )
            if not tile_geometry.intersects(raster_aoi, 1).getInfo():
                continue

            tile_id = f"r{row_index:+06d}_c{col_index:+06d}"
            tile_name = f"{CONFIG['direct_raster_prefix']}_{tile_id}"
            tile_path = raster_directory / f"{tile_name}.tif"
            tile_transform = [
                float(scale_m),
                0.0,
                float(tile_xmin),
                0.0,
                -float(scale_m),
                float(tile_ymax),
            ]
            download_url = raster_image.getDownloadURL(
                {
                    "bands": EXPECTED_BAND_NAMES,
                    "crs": ROUTE_L_CRS,
                    "crs_transform": tile_transform,
                    "dimensions": [tile_pixels, tile_pixels],
                    "format": "GEO_TIFF",
                }
            )
            _download_url(download_url, tile_path)
            if not _is_geotiff(tile_path):
                raise RuntimeError(f"Not a GeoTIFF: {tile_path}")
            direct_raster_manifest.append(
                {
                    "tile_id": tile_id,
                    "row_index": int(row_index),
                    "column_index": int(col_index),
                    "xmin": float(tile_xmin),
                    "ymin": float(tile_ymin),
                    "xmax": float(tile_xmax),
                    "ymax": float(tile_ymax),
                    "width_pixels": int(tile_pixels),
                    "height_pixels": int(tile_pixels),
                    "uncompressed_request_mb": round(
                        uncompressed_tile_mb, 3
                    ),
                    "file": str(tile_path.relative_to(direct_root)),
                    "bytes": int(tile_path.stat().st_size),
                    "sha256": _sha256(tile_path),
                }
            )
            print(
                f"Downloaded raster tile {len(direct_raster_manifest)}: "
                f"{tile_path.name}"
            )

    if not direct_raster_manifest:
        raise RuntimeError("No raster tiles intersected the configured AOI")
    pd.DataFrame(direct_raster_manifest).to_csv(
        direct_root / "raster_tile_manifest.csv", index=False
    )
    print(f"Downloaded {len(direct_raster_manifest)} AOI raster tiles")


# %% --------------------------------------------------------------------------
# 11. WRITE DIRECT MANIFEST/ZIP OR CHECK BATCH TASK STATES
# -----------------------------------------------------------------------------

if direct_root is not None:
    pd.DataFrame(direct_point_manifest).to_csv(
        direct_root / "point_download_manifest.csv", index=False
    )
    run_manifest = {
        "definition_version": "hierarchical_similarity_micro_pilot_v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "export_destination": EXPORT_DESTINATION,
        "point_shards": len(direct_point_manifest),
        "point_rows": sum(item["rows"] for item in direct_point_manifest),
        "raster_tiles": len(direct_raster_manifest),
        "raster_aoi_asset": CONFIG["direct_raster_aoi_asset"],
        "raster_tile_size_m": CONFIG["direct_raster_tile_size_m"],
        "route_l_crs": ROUTE_L_CRS,
        "route_l_transform_30m": ROUTE_L_TRANSFORM_30M,
        "missing_predictor_value": CONFIG["missing_predictor_value"],
        "qa_missing_value": 0,
        "colab_storage_warning": (
            "/content is temporary; download the ZIP before disconnecting"
        ),
    }
    with (direct_root / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(run_manifest, handle, indent=2)

    zip_parent = Path("/content") if IN_COLAB else Path("/tmp")
    zip_path = zip_parent / CONFIG["direct_zip_name"]
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(direct_root.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=path.relative_to(direct_root))
    print(f"Direct-download ZIP ready: {zip_path}")
    print(f"ZIP size: {zip_path.stat().st_size / (1024 ** 2):.1f} MB")
    if DOWNLOAD_DIRECT_ZIP_TO_BROWSER and IN_COLAB:
        print("Starting browser download. Keep this tab open until it finishes.")
        files.download(str(zip_path))

if export_tasks:
    for item in export_tasks:
        status = item["task"].status()
        print(
            item["database_code"],
            item["broad_ecosystem_group"],
            f"batch {item['batch']:02d}",
            f"expected rows {item['expected_rows']}",
            status.get("state"),
            status.get("error_message"),
        )
elif direct_root is None:
    print("No point-export tasks or direct downloads were started in this run")
