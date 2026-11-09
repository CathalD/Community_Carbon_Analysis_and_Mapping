"""Part 3: sample the existing AOI covariate stack at community cores.

Run Part 1 of Covariate_HBJB.ipynb first so ``covariate_stack`` and
``qa_stack`` already exist. This part reuses those AOI-clipped images and
starts a table export only; it does not rebuild the Canada-wide stack and it
does not submit any raster export.
"""


# %% --------------------------------------------------------------------------
# PART 3A. LOAD THE ONE-ROW-PER-COMMUNITY-CORE TABLE
# -----------------------------------------------------------------------------

import os

import ee
import pandas as pd

required_aoi_objects = [
    "CONFIG",
    "FOCAL_PROJECTION",
    "covariate_stack",
    "qa_stack",
    "expected_band_names",
    "export_region",
]
missing_aoi_objects = [
    name for name in required_aoi_objects if name not in globals()
]
if missing_aoi_objects:
    raise RuntimeError(
        "Run every Part 1 cell first. Missing AOI definitions: "
        + ", ".join(missing_aoi_objects)
    )
if len(expected_band_names) != 119:
    raise ValueError(
        f"Expected the 119-band AOI stack, found {len(expected_band_names)}."
    )

try:
    from google.colab import files
    IN_COLAB = True
except ImportError:
    IN_COLAB = False

# This is independent of the Part 1 raster and Part 2 external-core switches.
RUN_COMMUNITY_POINT_EXPORT = True
COMMUNITY_TABLE_LOCAL_PATH = "/content/community_cores_for_gee.csv"

if not os.path.exists(COMMUNITY_TABLE_LOCAL_PATH):
    if IN_COLAB:
        print("Upload community_cores_for_gee.csv")
        uploaded = files.upload()
        expected_name = "community_cores_for_gee.csv"
        if expected_name in uploaded:
            COMMUNITY_TABLE_LOCAL_PATH = f"/content/{expected_name}"
            if not os.path.exists(COMMUNITY_TABLE_LOCAL_PATH):
                COMMUNITY_TABLE_LOCAL_PATH = expected_name
        elif len(uploaded) == 1:
            uploaded_name = next(iter(uploaded))
            uploaded_path = f"/content/{uploaded_name}"
            COMMUNITY_TABLE_LOCAL_PATH = (
                uploaded_path
                if os.path.exists(uploaded_path)
                else uploaded_name
            )
        else:
            raise FileNotFoundError(
                "Could not identify community_cores_for_gee.csv in the "
                "uploaded files."
            )
    else:
        raise FileNotFoundError(
            f"{COMMUNITY_TABLE_LOCAL_PATH} not found. Upload the community "
            "CSV or change COMMUNITY_TABLE_LOCAL_PATH."
        )

required_community_cols = {
    "core_uid",
    "database_code",
    "plot_id",
    "latitude",
    "longitude",
    "country",
    "broad_ecosystem_group",
}
community_cores_df = pd.read_csv(
    COMMUNITY_TABLE_LOCAL_PATH, low_memory=False
)
missing_community_cols = required_community_cols - set(
    community_cores_df.columns
)
if missing_community_cols:
    raise ValueError(
        "Rebuild the community-core table before continuing. Missing "
        f"columns: {sorted(missing_community_cols)}"
    )

community_cores_df = community_cores_df[
    community_cores_df["country"].astype(str).str.casefold().eq("canada")
].copy()
community_cores_df = community_cores_df.dropna(
    subset=["latitude", "longitude"]
)
community_cores_df = community_cores_df[
    community_cores_df["latitude"].between(-90, 90)
    & community_cores_df["longitude"].between(-180, 180)
].drop_duplicates("core_uid")

print(
    f"Loaded {len(community_cores_df):,} unique community cores for Part 3"
)


# %% --------------------------------------------------------------------------
# PART 3B. BUILD AOI COMMUNITY POINTS AND SAMPLE THE PART 1 STACK
# -----------------------------------------------------------------------------


def _clean_community_property(value):
    """Convert pandas/numpy scalars to Earth Engine-safe Python values."""
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def _community_frame_to_features(frame):
    features = []
    for _, row in frame.iterrows():
        properties = {}
        for key, value in row.to_dict().items():
            clean_value = _clean_community_property(value)
            if clean_value is not None:
                properties[key] = clean_value
        point = ee.Geometry.Point(
            [float(row["longitude"]), float(row["latitude"])]
        )
        features.append(ee.Feature(point, properties))
    return ee.FeatureCollection(features)


# These eight reviewed community coordinates are already known to fall inside
# export_region, so retain them directly without a synchronous server check.
community_points = _community_frame_to_features(community_cores_df)
community_points_in_aoi = len(community_cores_df)

# Preserve all eight rows if an individual predictor is masked. The QA fields
# distinguish valid values from the -9999 predictor sentinel during analysis.
COMMUNITY_MISSING_PREDICTOR_VALUE = -9999
community_sampling_image = covariate_stack.unmask(
    COMMUNITY_MISSING_PREDICTOR_VALUE
).addBands(qa_stack.unmask(0))

community_samples = community_sampling_image.sampleRegions(
    collection=community_points,
    projection=FOCAL_PROJECTION,
    tileScale=8,
    geometries=True,
)

print(
    f"Community AOI sample staged: {community_points_in_aoi} cores, "
    "119 predictors + 5 QA fields"
)
if "Map" in globals():
    Map.addLayer(
        community_points,
        {"color": "orange"},
        "Part 3 community cores",
    )


# %% --------------------------------------------------------------------------
# PART 3C. EXPORT THE COMMUNITY-CORE AOI TABLE TO CLOUD STORAGE
# -----------------------------------------------------------------------------

COMMUNITY_GCS_BUCKET = "gee_exports_cd"
COMMUNITY_GCS_PREFIX = (
    f"{CONFIG['export_folder']}/CommunityCoreCovariatesAOI/"
    "RouteL_AOI_CommunityCoreCovariates"
)

community_export_task = None
if RUN_COMMUNITY_POINT_EXPORT:
    community_export_task = ee.batch.Export.table.toCloudStorage(
        collection=community_samples,
        description="RouteL_AOI_CommunityCoreCovariates_toGCS",
        bucket=COMMUNITY_GCS_BUCKET,
        fileNamePrefix=COMMUNITY_GCS_PREFIX,
        fileFormat="CSV",
    )
    community_export_task.start()
    print(
        f"Started community export {community_export_task.id}: "
        f"{COMMUNITY_GCS_PREFIX}.csv"
    )
else:
    print(
        "Community export skipped: set RUN_COMMUNITY_POINT_EXPORT = True "
        "when ready."
    )


# %% --------------------------------------------------------------------------
# PART 3D. CHECK THE COMMUNITY TASK WITHOUT RESTARTING IT
# -----------------------------------------------------------------------------

if community_export_task is None:
    print("No community-core task was submitted in this run.")
else:
    community_status = community_export_task.status()
    print(
        "community_cores",
        f"expected rows {community_points_in_aoi}",
        community_status.get("state"),
        community_status.get("error_message"),
    )
