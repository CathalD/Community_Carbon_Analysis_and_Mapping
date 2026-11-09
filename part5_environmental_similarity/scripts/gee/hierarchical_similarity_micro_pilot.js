/**
 * Hierarchical, multi-scale similarity micro-pilot (Earth Engine Code Editor).
 *
 * Purpose
 * -------
 * Build one compact, interpretable 30 m covariate stack for a deliberately
 * tiny pilot area, then export both the raster and the same covariates sampled
 * at configured PILOT/VALIDATION soil-core points inside that area. Every
 * predictor is prefixed by its intended
 * similarity block: local_, regional_, or broad_.  This is an isolated research
 * export; it does not replace or modify the authoritative Route L stack.
 *
 * Before running
 * --------------
 * 1. Change CONFIG.aoiAsset if needed.
 * 2. Upload a small pilot/validation point table as a FeatureCollection and set
 *    CONFIG.corePointsAsset. Retain a stable core_uid property. This script does
 *    NOT extract these new predictors for nationwide external donor cores.
 * 3. Inspect the default 2 km-wide AOI-centroid pilot.  To use a hand-drawn
 *    geometry, import it as `geometry`, set USE_DRAWN_GEOMETRY = true, and set
 *    DRAWN_PILOT_GEOMETRY = geometry below.
 * 4. Optionally configure river/coast/parent-material/permafrost assets.
 * 5. Run once with START_EXPORTS = false.  Inspect layers and printed bands,
 *    then switch it to true and run again.
 *
 * Nationwide donor note
 * ---------------------
 * External donors are intentionally out of scope here: the source collections
 * are spatially filtered for the tiny processing area, and the final stack is
 * clipped to the pilot. A nationwide donor-point extraction must rebuild the
 * identical definitions around donor points and export them in spatial shards
 * (adapt scripts/gee/canada_external_core_sampling.py). Do not treat this tiny
 * raster or its point CSV as a donor-covariate export.
 *
 * Deliberate exclusion: raw TWI is not included.  Earlier AOI testing found it
 * unstable in very low relief terrain; defensible hydrologic derivatives need
 * a conditioned DEM and an explicit flow-routing method outside this script.
 */


// -----------------------------------------------------------------------------
// 0. USER CONFIGURATION
// -----------------------------------------------------------------------------

var CONFIG = {
  aoiAsset: 'projects/north-star-project-470316/assets/aoi_fs',

  // Leave blank to preview/export the raster without sampling core points.
  // If set, only points inside pilotArea are sampled. This hook is for a few
  // pilot/validation points, not the nationwide external donor collection.
  corePointsAsset: '',

  // Tiny default: centroid of the configured AOI, buffered by this radius.
  // The exported rectangle is roughly 2 x radius on each side, clipped to AOI.
  pilotBufferM: 1000,

  startYear: 2019,
  endYear: 2023,
  s2CloudPercent: 40,
  scaleM: 30,
  localContextRadiusM: 120,
  regionalContextRadiusM: 1000,
  waterDistanceMaxM: 12000,

  // Optional vector hooks. Use line features for coastlines and rivers.
  // Empty strings keep the default run valid and simply omit these bands.
  coastlineAsset: '',
  riverAsset: '',
  optionalVectorDistanceMaxM: 250000,

  // Optional raster hooks. Set both asset and source-band name to include one.
  // These hooks assume class IDs. They are exported for exact-match or one-hot
  // handling downstream; never standardize the codes as continuous numbers.
  parentMaterialAsset: '',
  parentMaterialBand: '',
  permafrostAsset: '',
  permafrostBand: '',

  driveFolder: 'MultiScaleCarbonCovariates',
  exportPrefix: 'HierarchicalSimilarity_MicroPilot',
  missingPointValue: -9999
};

var START_EXPORTS = false;
var USE_DRAWN_GEOMETRY = false;
var DRAWN_PILOT_GEOMETRY = null; // e.g. set this to imported `geometry`.

if (CONFIG.endYear < CONFIG.startYear) {
  throw new Error('CONFIG.endYear must be >= CONFIG.startYear.');
}
if (USE_DRAWN_GEOMETRY && DRAWN_PILOT_GEOMETRY === null) {
  throw new Error('Import/draw a geometry before USE_DRAWN_GEOMETRY = true.');
}


// Exact 30 m Route L grid (custom Canada Lambert Conformal Conic / ESRI:102002).
// The transform preserves the same grid origin as the existing Route L exports.
var ROUTE_L_CRS =
  'PROJCS["Canada_Lambert_Conformal_Conic",' +
  'GEOGCS["GCS_North_American_1983",' +
  'DATUM["D_North_American_1983",' +
  'SPHEROID["GRS_1980",6378137,298.257222101]],' +
  'PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]],' +
  'PROJECTION["Lambert_Conformal_Conic_2SP"],' +
  'PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],' +
  'PARAMETER["Central_Meridian",-96],' +
  'PARAMETER["Standard_Parallel_1",50],' +
  'PARAMETER["Standard_Parallel_2",70],' +
  'PARAMETER["Latitude_Of_Origin",40],UNIT["Meter",1]]';

var ROUTE_L_TRANSFORM_30M = [
  30.0, 0.0, 462400.096580488,
  0.0, -30.0, 1881190.18705772
];
var OUTPUT_PROJECTION = ee.Projection(ROUTE_L_CRS, ROUTE_L_TRANSFORM_30M);


// -----------------------------------------------------------------------------
// 1. PILOT AND PROCESSING REGIONS
// -----------------------------------------------------------------------------

var aoiFeatures = ee.FeatureCollection(CONFIG.aoiAsset);
var aoi = aoiFeatures.geometry(1);
var defaultPilot = aoi.centroid(1)
  .buffer(CONFIG.pilotBufferM, 1)
  .bounds(1)
  .intersection(aoi, 1);
var pilotArea = USE_DRAWN_GEOMETRY ?
  ee.Geometry(DRAWN_PILOT_GEOMETRY).intersection(aoi, 1) : defaultPilot;

// Focal and distance calculations need context outside the final tiny export.
var contextBufferM = Math.max(
  CONFIG.regionalContextRadiusM + 60,
  CONFIG.waterDistanceMaxM + 60
);
var processingArea = pilotArea.buffer(contextBufferM, 1).intersection(aoi, 1);

Map.centerObject(pilotArea, 13);
Map.addLayer(
  ee.Image().byte().paint(aoiFeatures, 1, 2),
  {palette: ['d73027']},
  'Configured AOI',
  false
);
Map.addLayer(
  ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(pilotArea)]), 1, 3),
  {palette: ['ffff00']},
  'Micro-pilot export boundary'
);


// -----------------------------------------------------------------------------
// 2. HELPERS
// -----------------------------------------------------------------------------

function combineImages(images) {
  var result = ee.Image(images[0]);
  for (var i = 1; i < images.length; i++) {
    result = result.addBands(ee.Image(images[i]));
  }
  return result;
}

function maskSentinel2(image) {
  var scl = image.select('SCL');
  var clear = scl.eq(4).or(scl.eq(5)).or(scl.eq(6));
  return image
    .select(['B2', 'B3', 'B4', 'B8', 'B11'])
    .multiply(0.0001)
    .updateMask(clear)
    .copyProperties(image, ['system:time_start']);
}

function addSentinel2Indices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI');
  var evi = image.expression(
    '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
    {
      nir: image.select('B8'),
      red: image.select('B4'),
      blue: image.select('B2')
    }
  ).rename('EVI');
  return image.addBands([ndvi, ndmi, evi]);
}

function focalMean(image, radiusM, outputName) {
  return image.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.circle({radius: radiusM, units: 'meters'})
  }).rename(outputName);
}

function focalStdDev(image, radiusM, outputName) {
  return image.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel: ee.Kernel.circle({radius: radiusM, units: 'meters'})
  }).rename(outputName);
}

function addSentinel1Metrics(image) {
  // VV and VH are in dB. Their difference is the log-domain VV/VH ratio.
  var contrastDb = image.select('VV').subtract(image.select('VH'))
    .rename('VV_minus_VH_dB');
  var ratioLinear = ee.Image(10).pow(contrastDb.divide(10))
    .rename('VV_over_VH_linear');
  return image.select(['VV', 'VH'])
    .addBands([contrastDb, ratioLinear])
    .copyProperties(image, ['system:time_start']);
}

function optionalVectorDistance(assetId, outputName) {
  var features = ee.FeatureCollection(assetId)
    .filterBounds(processingArea.buffer(CONFIG.optionalVectorDistanceMaxM, 1));
  return features.distance({
    searchRadius: CONFIG.optionalVectorDistanceMaxM,
    maxError: 100
  }).unmask(CONFIG.optionalVectorDistanceMaxM).rename(outputName).toFloat();
}

function optionalRaster(assetId, sourceBand, outputName) {
  return ee.Image(assetId).select(sourceBand).rename(outputName).clip(processingArea);
}


// -----------------------------------------------------------------------------
// 3. LOCAL BLOCK: VEGETATION, WETNESS, RADAR AND FINE TERRAIN
// -----------------------------------------------------------------------------

var startDate = ee.Date.fromYMD(CONFIG.startYear, 1, 1);
var endDate = ee.Date.fromYMD(CONFIG.endYear + 1, 1, 1);

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(processingArea)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CONFIG.s2CloudPercent))
  .map(maskSentinel2)
  .map(addSentinel2Indices);

var s2Spring = s2.filter(ee.Filter.calendarRange(4, 6, 'month'));
var s2Summer = s2.filter(ee.Filter.calendarRange(6, 8, 'month'));

// Collection reductions otherwise inherit Earth Engine's generic WGS84
// default. Pin focal inputs to the shared 30 m carrier grid before computing
// neighbourhood summaries so raster and donor-point extractions are identical.
var ndmiSpringMedian = s2Spring.select('NDMI').median()
  .setDefaultProjection(OUTPUT_PROJECTION);
var ndviSummerMedian = s2Summer.select('NDVI').median()
  .setDefaultProjection(OUTPUT_PROJECTION);
var eviSummerMedian = s2Summer.select('EVI').median()
  .setDefaultProjection(OUTPUT_PROJECTION);

var localImages = [
  ndmiSpringMedian.rename('local_s2_ndmi_spring_median'),
  s2Spring.select('NDMI').reduce(ee.Reducer.stdDev())
    .rename('local_s2_ndmi_spring_temporal_sd'),
  ndviSummerMedian.rename('local_s2_ndvi_summer_median'),
  s2Summer.select('NDVI').reduce(ee.Reducer.stdDev())
    .rename('local_s2_ndvi_summer_temporal_sd'),
  eviSummerMedian.rename('local_s2_evi_summer_median'),
  s2Summer.select('EVI').reduce(ee.Reducer.stdDev())
    .rename('local_s2_evi_summer_temporal_sd'),
  focalMean(
    ndviSummerMedian, CONFIG.localContextRadiusM,
    'local_s2_ndvi_' + CONFIG.localContextRadiusM + 'm_mean'
  ),
  focalStdDev(
    ndviSummerMedian, CONFIG.localContextRadiusM,
    'local_s2_ndvi_' + CONFIG.localContextRadiusM + 'm_spatial_sd'
  )
];

var dynamicWorld = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterBounds(processingArea)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.calendarRange(6, 8, 'month'));

var dwBands = ['trees', 'shrub_and_scrub', 'flooded_vegetation', 'water'];
var dwNames = ['tree', 'shrub', 'floodedveg', 'water'];
for (var d = 0; d < dwBands.length; d++) {
  var dwMean = dynamicWorld.select(dwBands[d]).mean()
    .setDefaultProjection(OUTPUT_PROJECTION);
  localImages.push(dwMean.rename('local_dw_' + dwNames[d] + '_prob_mean'));
  localImages.push(
    dynamicWorld.select(dwBands[d]).reduce(ee.Reducer.stdDev())
      .rename('local_dw_' + dwNames[d] + '_prob_temporal_sd')
  );
  localImages.push(
    focalMean(
      dwMean, CONFIG.localContextRadiusM,
      'local_dw_' + dwNames[d] + '_' + CONFIG.localContextRadiusM + 'm_mean'
    )
  );
}

var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(processingArea)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .map(addSentinel1Metrics);

var s1MetricBands = ['VV', 'VH', 'VV_minus_VH_dB', 'VV_over_VH_linear'];
var s1MetricNames = ['vv_db', 'vh_db', 'vv_minus_vh_db', 'vv_over_vh_linear'];
for (var s = 0; s < s1MetricBands.length; s++) {
  localImages.push(
    sentinel1.select(s1MetricBands[s]).median()
      .rename('local_s1_' + s1MetricNames[s] + '_median')
  );
  localImages.push(
    sentinel1.select(s1MetricBands[s]).reduce(ee.Reducer.stdDev())
      .rename('local_s1_' + s1MetricNames[s] + '_temporal_sd')
  );
}

// CDEM descriptors only; no raw TWI.
var cdemCollection = ee.ImageCollection('NRCan/CDEM')
  .filterBounds(processingArea)
  .select('elevation');
var cdemFirst = ee.Image(cdemCollection.first());
var elevation = cdemCollection.mosaic()
  .setDefaultProjection(cdemFirst.projection())
  .resample('bilinear')
  .rename('elevation_m')
  .clip(processingArea);
var slope = ee.Terrain.slope(elevation);
var localElevationMean = focalMean(
  elevation, CONFIG.localContextRadiusM, 'local_elevation_context_mean'
);
localImages.push(slope.rename('local_slope_deg'));
localImages.push(
  elevation.subtract(localElevationMean)
    .rename('local_relative_elevation_' + CONFIG.localContextRadiusM + 'm')
);
localImages.push(
  focalStdDev(
    elevation, CONFIG.localContextRadiusM,
    'local_elevation_' + CONFIG.localContextRadiusM + 'm_sd'
  )
);

var localStack = combineImages(localImages).toFloat();


// -----------------------------------------------------------------------------
// 4. REGIONAL BLOCK: WATER REGIME AND MESOTOPOGRAPHY
// -----------------------------------------------------------------------------

var surfaceWater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
// JRC masks locations where surface water was never detected. For these three
// occurrence/frequency descriptors that state is a meaningful zero, not an
// unknown value. Filling it also prevents ordinary land pixels from losing the
// entire regional similarity block during complete-case sampling.
var regionalImages = [
  surfaceWater.select('occurrence')
    .unmask(0)
    .rename('regional_jrc_water_occurrence_pct'),
  surfaceWater.select('seasonality')
    .unmask(0)
    .rename('regional_jrc_water_seasonality_months'),
  surfaceWater.select('recurrence')
    .unmask(0)
    .rename('regional_jrc_water_recurrence_pct')
];

// A transparent proximity-to-observed-water proxy.  It is not equivalent to
// distance to a mapped river or coast; use the optional vector hooks for those.
var maxWaterDistancePixels = Math.ceil(
  CONFIG.waterDistanceMaxM / CONFIG.scaleM
);
var waterBinary30m = surfaceWater.select('max_extent').eq(1)
  .unmask(0).toByte().setDefaultProjection(OUTPUT_PROJECTION);
var waterDistance = waterBinary30m.fastDistanceTransform(
  maxWaterDistancePixels, 'pixels', 'squared_euclidean'
).sqrt().multiply(CONFIG.scaleM)
  .min(CONFIG.waterDistanceMaxM)
  .rename('regional_observed_water_distance_m');
regionalImages.push(waterDistance);

var regionalElevationMean = focalMean(
  elevation, CONFIG.regionalContextRadiusM, 'regional_elevation_context_mean'
);
regionalImages.push(elevation.rename('regional_elevation_m'));
regionalImages.push(
  elevation.subtract(regionalElevationMean)
    .rename('regional_relative_elevation_' + CONFIG.regionalContextRadiusM + 'm')
);
regionalImages.push(
  focalStdDev(
    elevation, CONFIG.regionalContextRadiusM,
    'regional_elevation_' + CONFIG.regionalContextRadiusM + 'm_sd'
  )
);

if (CONFIG.riverAsset !== '') {
  regionalImages.push(
    optionalVectorDistance(CONFIG.riverAsset, 'regional_river_distance_m')
  );
}
if (CONFIG.coastlineAsset !== '') {
  regionalImages.push(
    optionalVectorDistance(CONFIG.coastlineAsset, 'regional_coast_distance_m')
  );
}

var regionalStack = combineImages(regionalImages).toFloat();


// -----------------------------------------------------------------------------
// 5. BROAD BLOCK: MULTI-YEAR CLIMATE + OPTIONAL SUBSTRATE/PERMAFROST
// -----------------------------------------------------------------------------

var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterBounds(processingArea)
  .filterDate(startDate, endDate);
var years = ee.List.sequence(CONFIG.startYear, CONFIG.endYear);

var annualPrecipitation = ee.ImageCollection.fromImages(years.map(function(y) {
  y = ee.Number(y);
  var yearStart = ee.Date.fromYMD(y, 1, 1);
  return era5.filterDate(yearStart, yearStart.advance(1, 'year'))
    .select('total_precipitation_sum')
    .sum()
    .multiply(1000)
    .rename('annual_precip_mm')
    .set('year', y);
}));

var annualFreezingMonths = ee.ImageCollection.fromImages(years.map(function(y) {
  y = ee.Number(y);
  var yearStart = ee.Date.fromYMD(y, 1, 1);
  return era5.filterDate(yearStart, yearStart.advance(1, 'year'))
    .select('temperature_2m')
    .map(function(image) { return image.lt(273.15); })
    .sum()
    .rename('freezing_months')
    .set('year', y);
}));

var broadImages = [
  era5.select('temperature_2m').mean().subtract(273.15)
    .rename('broad_era5_mean_temperature_c'),
  era5.select('temperature_2m').reduce(ee.Reducer.stdDev())
    .rename('broad_era5_monthly_temperature_sd_c'),
  annualPrecipitation.mean()
    .rename('broad_era5_mean_annual_precip_mm'),
  annualPrecipitation.reduce(ee.Reducer.stdDev())
    .rename('broad_era5_annual_precip_sd_mm'),
  annualFreezingMonths.mean()
    .rename('broad_era5_mean_freezing_months')
];

if (CONFIG.parentMaterialAsset !== '' && CONFIG.parentMaterialBand !== '') {
  broadImages.push(optionalRaster(
    CONFIG.parentMaterialAsset,
    CONFIG.parentMaterialBand,
    'broad_parent_material_categorical_code'
  ));
}
if (CONFIG.permafrostAsset !== '' && CONFIG.permafrostBand !== '') {
  broadImages.push(optionalRaster(
    CONFIG.permafrostAsset,
    CONFIG.permafrostBand,
    'broad_permafrost_or_groundice_categorical_code'
  ));
}

// IMPORTANT: the two optional *_categorical_code bands, when enabled, must be
// represented as one-hot indicators or an exact-match similarity term in the
// downstream model. Raw class-number differences have no metric meaning.

var broadStack = combineImages(broadImages).toFloat();


// -----------------------------------------------------------------------------
// 6. FINAL STACK, QA AND MAP REVIEW
// -----------------------------------------------------------------------------

// Cross-language contract shared with hierarchical_similarity_donor_sampling.py.
// Keep this order stable: the R workflow matches by name, and the donor sampler
// performs an exact server-side name/order check before it can submit tasks.
var DEFAULT_BAND_NAMES = [
  'local_s2_ndmi_spring_median',
  'local_s2_ndmi_spring_temporal_sd',
  'local_s2_ndvi_summer_median',
  'local_s2_ndvi_summer_temporal_sd',
  'local_s2_evi_summer_median',
  'local_s2_evi_summer_temporal_sd',
  'local_s2_ndvi_120m_mean',
  'local_s2_ndvi_120m_spatial_sd',
  'local_dw_tree_prob_mean',
  'local_dw_tree_prob_temporal_sd',
  'local_dw_tree_120m_mean',
  'local_dw_shrub_prob_mean',
  'local_dw_shrub_prob_temporal_sd',
  'local_dw_shrub_120m_mean',
  'local_dw_floodedveg_prob_mean',
  'local_dw_floodedveg_prob_temporal_sd',
  'local_dw_floodedveg_120m_mean',
  'local_dw_water_prob_mean',
  'local_dw_water_prob_temporal_sd',
  'local_dw_water_120m_mean',
  'local_s1_vv_db_median',
  'local_s1_vv_db_temporal_sd',
  'local_s1_vh_db_median',
  'local_s1_vh_db_temporal_sd',
  'local_s1_vv_minus_vh_db_median',
  'local_s1_vv_minus_vh_db_temporal_sd',
  'local_s1_vv_over_vh_linear_median',
  'local_s1_vv_over_vh_linear_temporal_sd',
  'local_slope_deg',
  'local_relative_elevation_120m',
  'local_elevation_120m_sd',
  'regional_jrc_water_occurrence_pct',
  'regional_jrc_water_seasonality_months',
  'regional_jrc_water_recurrence_pct',
  'regional_observed_water_distance_m',
  'regional_elevation_m',
  'regional_relative_elevation_1000m',
  'regional_elevation_1000m_sd',
  'broad_era5_mean_temperature_c',
  'broad_era5_monthly_temperature_sd_c',
  'broad_era5_mean_annual_precip_mm',
  'broad_era5_annual_precip_sd_mm',
  'broad_era5_mean_freezing_months',
  'qa_s2_spring_observations',
  'qa_s2_summer_observations',
  'qa_dynamic_world_observations',
  'qa_sentinel1_observations'
];
if (DEFAULT_BAND_NAMES.length !== 47) {
  throw new Error('The default cross-language band contract must contain 47 bands.');
}

var qaStack = combineImages([
  s2Spring.select('NDMI').count().rename('qa_s2_spring_observations'),
  s2Summer.select('NDVI').count().rename('qa_s2_summer_observations'),
  dynamicWorld.select('trees').count().rename('qa_dynamic_world_observations'),
  sentinel1.select('VV').count().rename('qa_sentinel1_observations')
]).unmask(0).toFloat();

var selectedBandNames = DEFAULT_BAND_NAMES.slice();
if (CONFIG.riverAsset !== '') {
  selectedBandNames.push('regional_river_distance_m');
}
if (CONFIG.coastlineAsset !== '') {
  selectedBandNames.push('regional_coast_distance_m');
}
if (CONFIG.parentMaterialAsset !== '' && CONFIG.parentMaterialBand !== '') {
  selectedBandNames.push('broad_parent_material_categorical_code');
}
if (CONFIG.permafrostAsset !== '' && CONFIG.permafrostBand !== '') {
  selectedBandNames.push('broad_permafrost_or_groundice_categorical_code');
}

var covariateStack = localStack
  .addBands(regionalStack)
  .addBands(broadStack)
  .addBands(qaStack)
  .select(selectedBandNames)
  .clip(pilotArea)
  .set({
    pilot_name: CONFIG.exportPrefix,
    purpose: 'hierarchical multi-scale similarity micro-pilot',
    scale_families: 'local|regional|broad',
    start_year: CONFIG.startYear,
    end_year: CONFIG.endYear,
    output_scale_m: CONFIG.scaleM,
    point_sampling_scope: 'pilot validation points inside raster only',
    raw_twi_included: false,
    categorical_band_rule: 'one-hot or exact-match; never numeric distance',
    route_l_authoritative_output: false
  });

print('Pilot area (ha)', pilotArea.area(1).divide(10000));
print('Covariate band count', covariateStack.bandNames().size());
print('Covariate bands', covariateStack.bandNames());
print('Raster metadata', covariateStack.toDictionary());

Map.addLayer(
  covariateStack.select('local_s2_ndvi_summer_median'),
  {min: -0.1, max: 0.9, palette: ['8c510a', 'f6e8c3', '01665e']},
  'Local: summer NDVI'
);
Map.addLayer(
  covariateStack.select('regional_observed_water_distance_m'),
  {min: 0, max: CONFIG.waterDistanceMaxM,
    palette: ['08306b', '6baed6', 'f7fbff']},
  'Regional: distance to observed water',
  false
);
Map.addLayer(
  covariateStack.select('broad_era5_mean_annual_precip_mm'),
  {min: 300, max: 1000, palette: ['f7fcf5', '74c476', '00441b']},
  'Broad: mean annual precipitation',
  false
);


// -----------------------------------------------------------------------------
// 7. OPTIONAL PILOT/VALIDATION-POINT SAMPLING
// -----------------------------------------------------------------------------

var hasCoreAsset = CONFIG.corePointsAsset !== '';
var pilotCores = ee.FeatureCollection([]);
var coreSamples = ee.FeatureCollection([]);

if (hasCoreAsset) {
  // The explicit spatial filter is a safety boundary. Nationwide external
  // donor extraction must use a separate, sharded workflow built from the same
  // band definitions; this pilot-clipped image cannot provide donor covariates.
  pilotCores = ee.FeatureCollection(CONFIG.corePointsAsset).filterBounds(pilotArea);
  // Unmask only for the point table so a missing predictor does not silently
  // remove an entire core.  QA fields identify weak observation support.
  coreSamples = covariateStack.unmask(CONFIG.missingPointValue).sampleRegions({
    collection: pilotCores,
    projection: OUTPUT_PROJECTION,
    tileScale: 4,
    geometries: true
  });
  print('Core points inside pilot', pilotCores.size());
  print('Point scope: pilot/validation only; nationwide donors are not sampled.');
  print('Sampled core preview', coreSamples.limit(10));
  Map.addLayer(pilotCores, {color: 'ff7f00'}, 'Configured pilot cores');
} else {
  print('Core sampling skipped: set CONFIG.corePointsAsset to a point asset.');
}


// -----------------------------------------------------------------------------
// 8. EXPORTS (TASKS ARE CREATED ONLY WHEN START_EXPORTS = true)
// -----------------------------------------------------------------------------

if (START_EXPORTS) {
  Export.image.toDrive({
    image: covariateStack.toFloat(),
    description: CONFIG.exportPrefix + '_Raster_30m',
    folder: CONFIG.driveFolder,
    fileNamePrefix: CONFIG.exportPrefix + '_Raster_30m',
    region: pilotArea,
    crs: ROUTE_L_CRS,
    crsTransform: ROUTE_L_TRANSFORM_30M,
    maxPixels: 1e10,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });

  if (hasCoreAsset) {
    Export.table.toDrive({
      collection: coreSamples,
      description: CONFIG.exportPrefix + '_CoreCovariates',
      folder: CONFIG.driveFolder,
      fileNamePrefix: CONFIG.exportPrefix + '_CoreCovariates',
      fileFormat: 'CSV'
    });
  } else {
    print('Raster task created; point task omitted because corePointsAsset is blank.');
  }
} else {
  print('No export tasks created. Inspect first, then set START_EXPORTS = true.');
}
