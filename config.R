# Settings for the workshop adaptation. No Route L models or outputs are sourced.
fs_config <- list(
  analysis_crs = "EPSG:3978", analysis_resolution_m = 1000,
  range_m = 2000, observation_log_sd = 0,
  sothe_uncertainty = "S100 full 90% interval width / (2*qnorm(.95)); absolute S30 sensitivity proxy",
  li_depth_observation_sd_cm = 5, # explicit working assumption; sensitivity-tested
  li_contact_mode = "all_eight_layer_evidence",
  li_organic_soc_threshold_g_kg = 80,
  li_suite_update_scale = "original_scale_gaussian", # retains genuine zero peat observations
  source_root = "..",
  sothe_root = "/Users/cathaldoherty/Desktop/Carbon_Resources/McMaster_WWFCanada_soil_carbon_250m_kg-m2_version3.0",
  covariate_root = "../data/raw/gee_imports/RouteL_Covariates_S2_2024_2025_AE_2024",
  # Pre-specified small models: do not fit 119 coefficients to eight cores.
  covariate_models = list(Mean_only = character(), Elevation = "elevation_m",
    Slope = "slope_deg", Relative_elevation = "relative_elevation_250m",
    Elevation_variability = "elevation_sd_250m", NDVI = "NDVI_summer",
    EVI = "EVI_summer", Trees = "DW_tree_prob", Shrubs = "DW_shrub_prob",
    Flooded_vegetation = "DW_floodedveg_prob", Wetness = "NDMI_spring",
    Combined = c("elevation_m", "NDVI_summer"),
    AlphaEarth_00 = "AE_00", AlphaEarth_01 = "AE_01"))
