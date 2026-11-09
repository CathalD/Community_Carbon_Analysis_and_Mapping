"""Build a multi-scale environmental covariate stack in Earth Engine.

Run this file in Google Colab. It acquires and visualizes covariates only: it
does not read community-core coordinates, fit a carbon model, calculate
environmental similarity, or modify Route L outputs.

The default run uses a central 10 km test area. Review the map, collection
counts, masks and test exports before setting CONFIG["test_only"] to False.
"""

# In a fresh Colab runtime, run once:
# !pip install geemap --quiet

import ee
import geemap


# -----------------------------------------------------------------------------
# 0. AUTHENTICATION AND CONFIGURATION
# -----------------------------------------------------------------------------

GEE_PROJECT = "north-star-project-470316"

ee.Authenticate()
ee.Initialize(project=GEE_PROJECT)

Map = geemap.Map()
Map.add_basemap("SATELLITE")

CONFIG = {
    "aoi_asset": "projects/north-star-project-470316/assets/aoi_fs",
    # Li products are visual references only; they are not predictors below.
    "li_carbon_stock_asset": (
        "projects/ee-cathalpdoherty2/assets/McMasterCarbon30mkgm2version1"
    ),
    "li_carbon_uncertainty_asset": (
        "projects/ee-cathalpdoherty2/assets/McMasterUncertaintyCarbon30mkgm2"
    ),
    "li_peat_depth_asset": (
        "projects/northstarlabs/assets/McMaster_PeatDepth_30m"
    ),
    "li_peat_depth_uncertainty_asset": (
        "projects/northstarlabs/assets/McMaster_Uncertainty_PeatDepth_30m"
    ),
    # These pooled years are provisional until field dates are confirmed.
    "years": [2019, 2020, 2021, 2022, 2023],
    "spring_window": ((4, 15), (6, 15)),
    "summer_window": ((6, 16), (8, 31)),
    "s2_cloud_threshold": 40,
    "min_s2_observations": 5,
    "alphaearth_year": 2021,
    "focal_radii_m": [30, 100, 250],
    "include_terrain_context": True,
    "terrain_context_radius_m": 250,
    # Keep this True until every small export has been downloaded and checked.
    "test_only": True,
    "test_half_width_m": 5000,
    "export_folder": "MultiScaleCarbonCovariates",
    "export_prefix": "RouteL_Covariates",
    # 8192 is a multiple of the default 256-pixel shard size.
    "file_dimensions": 8192,
}

if CONFIG["years"] != list(range(CONFIG["years"][0], CONFIG["years"][-1] + 1)):
    raise ValueError("CONFIG['years'] must be a contiguous ascending range.")
if CONFIG["alphaearth_year"] is None:
    raise ValueError("Choose one explicit AlphaEarth year for this comparison.")


# Exact Route L projected grid. The CRS is the Li product's custom Canada
# Lambert Conformal Conic (ESRI:102002), not EPSG:3347. The 10 m transform uses
# the same origin as the production 30 m grid, so each 3 x 3 block nests exactly
# inside a Route L cell. Export parameters, not Image.reproject(), set the grid.
ROUTE_L_GRID = {
    "crs": (
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
    ),
    "transform_10m": [
        10.0,
        0.0,
        462400.096580488,
        0.0,
        -10.0,
        1881190.18705772,
    ],
}

FOCAL_PROJECTION = ee.Projection(
    ROUTE_L_GRID["crs"], ROUTE_L_GRID["transform_10m"]
)


# -----------------------------------------------------------------------------
# 1. AOI, TEST REGION AND LI REFERENCE LAYERS
# -----------------------------------------------------------------------------

aoi_features = ee.FeatureCollection(CONFIG["aoi_asset"])
aoi = aoi_features.geometry(maxError=1)

if CONFIG["test_only"]:
    export_region = (
        aoi.centroid(maxError=1)
        .buffer(CONFIG["test_half_width_m"], maxError=1)
        .bounds(maxError=1)
        .intersection(aoi, maxError=1)
    )
else:
    export_region = aoi

max_focal_radius = max(CONFIG["focal_radii_m"])
processing_region = export_region.buffer(max_focal_radius + 60, maxError=1)

Map.centerObject(aoi, 9)
Map.addLayer(
    ee.Image().paint(aoi_features, 1, 2), {"palette": ["red"]}, "AOI boundary"
)
Map.addLayer(
    ee.Image().paint(export_region, 1, 2),
    {"palette": ["yellow"]},
    "Current export region",
)

li_carbon_stock = ee.Image(CONFIG["li_carbon_stock_asset"]).clip(aoi)
li_carbon_uncertainty = ee.Image(
    CONFIG["li_carbon_uncertainty_asset"]
).clip(aoi)
li_peat_depth = ee.Image(CONFIG["li_peat_depth_asset"]).clip(aoi)
li_peat_depth_uncertainty = ee.Image(
    CONFIG["li_peat_depth_uncertainty_asset"]
).clip(aoi)

Map.addLayer(
    li_carbon_stock,
    {"min": 0, "max": 30, "palette": ["f7fbff", "6baed6", "08306b"]},
    "Li carbon stock (reference only)",
    False,
)
Map.addLayer(
    li_carbon_uncertainty,
    {"min": 0, "max": 10, "palette": ["ffffcc", "fd8d3c", "bd0026"]},
    "Li carbon spread (reference only)",
    False,
)
Map.addLayer(
    li_peat_depth,
    {"min": 0, "max": 300, "palette": ["f7fcf5", "74c476", "00441b"]},
    "Li peat depth (reference only)",
    False,
)
Map.addLayer(
    li_peat_depth_uncertainty,
    {"min": 0, "max": 100, "palette": ["ffffcc", "fd8d3c", "bd0026"]},
    "Li depth spread (reference only)",
    False,
)


# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------

def fixed_calendar_window(start_month, start_day, end_month, end_day):
    """Build one inclusive month/day filter without per-year merges.

    Calendar fields preserve April 15 as April 15 in leap years, unlike a
    fixed day-of-year range.
    """
    if (end_month, end_day) < (start_month, start_day):
        raise ValueError("Seasonal windows may not cross a calendar year.")

    pieces = []
    if start_month == end_month:
        return ee.Filter.And(
            ee.Filter.calendarRange(start_month, start_month, "month"),
            ee.Filter.calendarRange(start_day, end_day, "day_of_month"),
        )

    pieces.append(
        ee.Filter.And(
            ee.Filter.calendarRange(start_month, start_month, "month"),
            ee.Filter.calendarRange(start_day, 31, "day_of_month"),
        )
    )
    if end_month > start_month + 1:
        pieces.append(
            ee.Filter.calendarRange(start_month + 1, end_month - 1, "month")
        )
    pieces.append(
        ee.Filter.And(
            ee.Filter.calendarRange(end_month, end_month, "month"),
            ee.Filter.calendarRange(1, end_day, "day_of_month"),
        )
    )
    return ee.Filter.Or(*pieces)


def mask_s2_clouds(image):
    """Mask cloud, shadow, snow and defective pixels with Sentinel-2 SCL."""
    scl = image.select("SCL")
    clear = scl.eq(4).Or(scl.eq(5)).Or(scl.eq(6))
    reflectance = (
        image.select(["B2", "B3", "B4", "B8", "B11", "B12"])
        .multiply(0.0001)
        .updateMask(clear)
    )
    return reflectance.copyProperties(image, ["system:time_start"])


def add_spring_wetness_indices(image):
    """Add NIR-SWIR wetness and Sentinel-2 tasseled-cap wetness."""
    ndmi = image.normalizedDifference(["B8", "B11"]).rename("NDMI_spring")
    tcw = image.expression(
        (
            "0.1509*B2 + 0.1973*B3 + 0.3279*B4 + 0.3406*B8 "
            "- 0.7112*B11 - 0.4572*B12"
        ),
        {
            "B2": image.select("B2"),
            "B3": image.select("B3"),
            "B4": image.select("B4"),
            "B8": image.select("B8"),
            "B11": image.select("B11"),
            "B12": image.select("B12"),
        },
    ).rename("TCW_spring")
    return image.addBands([ndmi, tcw])


def add_summer_greenness_indices(image):
    """Add summer NDVI and EVI."""
    ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI_summer")
    evi = image.expression(
        "2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1))",
        {
            "NIR": image.select("B8"),
            "RED": image.select("B4"),
            "BLUE": image.select("B2"),
        },
    ).rename("EVI_summer")
    return image.addBands([ndvi, evi])


def add_focal_stats(image, band_names, radii_m):
    """Append focal means and SDs for explicit circular radii in metres."""
    base = image.select(band_names)
    focal_images = []
    for radius in radii_m:
        kernel = ee.Kernel.circle(radius=radius, units="meters")
        focal_images.extend(
            [
                base.reduceNeighborhood(ee.Reducer.mean(), kernel).rename(
                    [f"{band}_focalMean_{radius}m" for band in band_names]
                ),
                base.reduceNeighborhood(ee.Reducer.stdDev(), kernel).rename(
                    [f"{band}_focalStd_{radius}m" for band in band_names]
                ),
            ]
        )
    return ee.Image.cat(focal_images)


# -----------------------------------------------------------------------------
# 3. SEASONAL SENTINEL-2 COMPOSITES AND QA
# -----------------------------------------------------------------------------

year_filter = ee.Filter.calendarRange(
    CONFIG["years"][0], CONFIG["years"][-1], "year"
)
spring_filter = fixed_calendar_window(
    *CONFIG["spring_window"][0], *CONFIG["spring_window"][1]
)
summer_filter = fixed_calendar_window(
    *CONFIG["summer_window"][0], *CONFIG["summer_window"][1]
)

s2_collection = (
    ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(processing_region)
    .filter(year_filter)
    .filter(
        ee.Filter.lt(
            "CLOUDY_PIXEL_PERCENTAGE", CONFIG["s2_cloud_threshold"]
        )
    )
)

s2_spring = (
    s2_collection.filter(spring_filter)
    .map(mask_s2_clouds)
    .map(add_spring_wetness_indices)
)
s2_summer = (
    s2_collection.filter(summer_filter)
    .map(mask_s2_clouds)
    .map(add_summer_greenness_indices)
)

spring_count = s2_spring.select("NDMI_spring").count().rename("S2_spring_count")
summer_count = s2_summer.select("NDVI_summer").count().rename("S2_summer_count")

spring_composite = (
    s2_spring.select(["NDMI_spring", "TCW_spring"])
    .median()
    .setDefaultProjection(FOCAL_PROJECTION)
    .updateMask(spring_count.gte(CONFIG["min_s2_observations"]))
    .clip(processing_region)
)
summer_composite = (
    s2_summer.select(["NDVI_summer", "EVI_summer"])
    .median()
    .setDefaultProjection(FOCAL_PROJECTION)
    .updateMask(summer_count.gte(CONFIG["min_s2_observations"]))
    .clip(processing_region)
)

Map.addLayer(
    spring_composite.select("TCW_spring"),
    {"min": -0.3, "max": 0.3, "palette": ["8c510a", "f6e8c3", "01665e"]},
    "Spring wetness (TCW)",
)
Map.addLayer(
    summer_composite.select("NDVI_summer"),
    {"min": -0.1, "max": 0.9, "palette": ["d73027", "fee08b", "1a9850"]},
    "Summer greenness (NDVI)",
)


# -----------------------------------------------------------------------------
# 4. DYNAMIC WORLD SOFT CLASS PROBABILITIES
# -----------------------------------------------------------------------------

# Use the same summer window as the greenness composite. These are model
# probabilities, not measured cover fractions, peat depth or carbon confidence.
dw_collection = (
    ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
    .filterBounds(processing_region)
    .filter(year_filter)
    .filter(summer_filter)
)
dw_source_bands = ["trees", "shrub_and_scrub", "flooded_vegetation"]
dw_output_bands = ["DW_tree_prob", "DW_shrub_prob", "DW_floodedveg_prob"]
dw_count = dw_collection.select("trees").count().rename("DW_summer_count")
dw_probs = (
    dw_collection.select(dw_source_bands)
    .mean()
    .rename(dw_output_bands)
    .setDefaultProjection(FOCAL_PROJECTION)
    .clip(processing_region)
)

Map.addLayer(
    dw_probs.select("DW_floodedveg_prob"),
    {"min": 0, "max": 1, "palette": ["ffffff", "4292c6", "08306b"]},
    "Dynamic World: flooded vegetation probability",
)


# -----------------------------------------------------------------------------
# 5. OPTIONAL TERRAIN CONTEXT
# -----------------------------------------------------------------------------

# The earlier focal-count/slope expression was not upslope contributing area
# and therefore was not a defensible TWI. This block exports transparent terrain
# descriptors instead. A true TWI would require a conditioned DEM plus a named
# flow-direction/accumulation method and should be a separate hypothesis.
terrain_stack = None
terrain_band_names = []
if CONFIG["include_terrain_context"]:
    dem_collection = (
        ee.ImageCollection("NRCan/CDEM")
        .filterBounds(processing_region)
        .select("elevation")
    )
    dem_first = ee.Image(dem_collection.first())
    dem = (
        dem_collection.mosaic()
        .setDefaultProjection(dem_first.projection())
        .rename("elevation_m")
        .clip(processing_region)
    )
    slope = ee.Terrain.slope(dem).rename("slope_deg")
    aspect_rad = ee.Terrain.aspect(dem).multiply(3.141592653589793 / 180)
    eastness = aspect_rad.sin().rename("eastness")
    northness = aspect_rad.cos().rename("northness")
    terrain_kernel = ee.Kernel.circle(
        radius=CONFIG["terrain_context_radius_m"], units="meters"
    )
    local_mean = dem.reduceNeighborhood(ee.Reducer.mean(), terrain_kernel)
    relative_elevation = dem.subtract(local_mean).rename(
        f"relative_elevation_{CONFIG['terrain_context_radius_m']}m"
    )
    elevation_sd = dem.reduceNeighborhood(
        ee.Reducer.stdDev(), terrain_kernel
    ).rename(f"elevation_sd_{CONFIG['terrain_context_radius_m']}m")
    terrain_stack = ee.Image.cat(
        [dem, slope, eastness, northness, relative_elevation, elevation_sd]
    )
    terrain_band_names = [
        "elevation_m",
        "slope_deg",
        "eastness",
        "northness",
        f"relative_elevation_{CONFIG['terrain_context_radius_m']}m",
        f"elevation_sd_{CONFIG['terrain_context_radius_m']}m",
    ]


# -----------------------------------------------------------------------------
# 6. ALPHAEARTH EMBEDDINGS
# -----------------------------------------------------------------------------

alpha_source_bands = [f"A{axis:02d}" for axis in range(64)]
alpha_output_bands = [f"AE_{axis:02d}" for axis in range(64)]
alpha_start = f"{CONFIG['alphaearth_year']}-01-01"
alpha_end = f"{CONFIG['alphaearth_year'] + 1}-01-01"

# The collection is spatially tiled. Mosaic one explicit annual layer; do not
# take a temporal median of axes whose joint vector is the meaningful object.
alphaearth_img = (
    ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")
    .filterBounds(processing_region)
    .filterDate(alpha_start, alpha_end)
    .select(alpha_source_bands)
    .mosaic()
    .rename(alpha_output_bands)
    .clip(processing_region)
)

Map.addLayer(
    alphaearth_img.select(["AE_01", "AE_16", "AE_09"]),
    {"min": -0.3, "max": 0.3},
    "AlphaEarth embedding preview",
    False,
)


# -----------------------------------------------------------------------------
# 7. MULTI-SCALE SUMMARIES AND FINAL STACKS
# -----------------------------------------------------------------------------

focal_source_bands = [
    "NDMI_spring",
    "TCW_spring",
    "NDVI_summer",
    "EVI_summer",
    "DW_tree_prob",
    "DW_shrub_prob",
    "DW_floodedveg_prob",
]
focal_source = spring_composite.addBands(summer_composite).addBands(dw_probs)
focal_stats_stack = add_focal_stats(
    focal_source, focal_source_bands, CONFIG["focal_radii_m"]
)

focal_band_names = []
for radius in CONFIG["focal_radii_m"]:
    focal_band_names.extend(
        [f"{band}_focalMean_{radius}m" for band in focal_source_bands]
    )
    focal_band_names.extend(
        [f"{band}_focalStd_{radius}m" for band in focal_source_bands]
    )

surface_context = (
    focal_source.addBands(focal_stats_stack)
    .select(focal_source_bands + focal_band_names)
    .clip(export_region)
)
embeddings = alphaearth_img.select(alpha_output_bands).clip(export_region)

stack_components = [surface_context, embeddings]
if terrain_stack is not None:
    terrain_stack = terrain_stack.select(terrain_band_names).clip(export_region)
    stack_components.append(terrain_stack)
covariate_stack = ee.Image.cat(stack_components).clip(export_region)

qa_stack = (
    ee.Image.cat(
        [
            spring_count,
            summer_count,
            dw_count,
            surface_context.mask().reduce(ee.Reducer.min()).rename(
                "surface_context_complete"
            ),
            embeddings.mask().reduce(ee.Reducer.min()).rename(
                "embeddings_available"
            ),
        ]
    )
    .unmask(0)
    .clip(export_region)
)

expected_band_names = (
    focal_source_bands + focal_band_names + alpha_output_bands + terrain_band_names
)
print("Expected covariate bands:", len(expected_band_names))
print(expected_band_names)
print(
    "Test export enabled."
    if CONFIG["test_only"]
    else "FULL AOI selected: confirm test outputs before starting tasks."
)

Map.addLayer(
    focal_stats_stack.select("NDVI_summer_focalMean_250m"),
    {"min": -0.1, "max": 0.9, "palette": ["d73027", "fee08b", "1a9850"]},
    "NDVI summer: 250 m focal mean",
    False,
)
Map.addLayer(
    qa_stack.select("S2_summer_count"),
    {"min": 0, "max": 30},
    "Sentinel-2 summer observation count",
    False,
)


# -----------------------------------------------------------------------------
# 8. MODULAR ASYNCHRONOUS EXPORTS
# -----------------------------------------------------------------------------

export_groups = {
    "surface_context": surface_context,
    "embeddings": embeddings,
    "QA": qa_stack,
}
if terrain_stack is not None:
    export_groups["terrain"] = terrain_stack


def export_kwargs(image, description, region):
    """Common export settings with exact Route L-aligned 10 m pixels."""
    return {
        "image": image.toFloat(),
        "description": description,
        "region": region,
        "crs": ROUTE_L_GRID["crs"],
        "crsTransform": ROUTE_L_GRID["transform_10m"],
        "maxPixels": 1e13,
        "fileFormat": "GeoTIFF",
        "fileDimensions": CONFIG["file_dimensions"],
        "skipEmptyTiles": True,
        "formatOptions": {"cloudOptimized": True},
    }


def export_groups_to_drive(groups, region, config):
    """Start one Drive task per logical group and return the tasks."""
    extent = "TEST" if config["test_only"] else "AOI"
    tasks = []
    for group_name, image in groups.items():
        prefix = f"{config['export_prefix']}_{extent}_{group_name}"
        kwargs = export_kwargs(image, f"{prefix}_toDrive", region)
        kwargs.update(
            {"folder": config["export_folder"], "fileNamePrefix": prefix}
        )
        task = ee.batch.Export.image.toDrive(**kwargs)
        task.start()
        tasks.append(task)
        print(f"Started Drive task {task.id}: {prefix}")
    return tasks


def export_groups_to_gcs(groups, region, config, bucket):
    """Start one Cloud Storage task per logical group and return the tasks."""
    extent = "TEST" if config["test_only"] else "AOI"
    tasks = []
    for group_name, image in groups.items():
        prefix = f"{config['export_prefix']}_{extent}_{group_name}"
        kwargs = export_kwargs(image, f"{prefix}_toGCS", region)
        kwargs.update(
            {
                "bucket": bucket,
                "fileNamePrefix": f"{config['export_folder']}/{prefix}",
            }
        )
        task = ee.batch.Export.image.toCloudStorage(**kwargs)
        task.start()
        tasks.append(task)
        print(f"Started GCS task {task.id}: {prefix}")
    return tasks


# Start only one destination after inspecting the map and QA layers:
# drive_tasks = export_groups_to_drive(export_groups, export_region, CONFIG)
# gcs_tasks = export_groups_to_gcs(
#     export_groups, export_region, CONFIG, bucket="your-gcs-bucket"
# )

Map
