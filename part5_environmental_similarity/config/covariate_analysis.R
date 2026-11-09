# Locked protocol for the first Route L covariate experiment.
#
# This configuration records scientific decisions. It does not modify or
# overwrite the production Route L configuration in config/route_l.R.
list(
  plan_id = "route_l_covariate_experiment_v1",
  status = "locked",
  locked_on = "2026-09-21",
  clarified_on = "2026-09-23",

  baseline = list(
    role = "authoritative_comparator__never_overwrite_during_candidate_tests",
    map_grid_m = 30,
    geographic_weight = "exp(-distance_m / bandwidth_m)",
    primary_bandwidth_m = 2000,
    uncertainty_label = "conditional_SD_proxy__not_a_calibrated_posterior",
    mineral_model = "retain_current_simple_distance_weighted_model"
  ),

  local_similarity_experiment = list(
    geographic_bandwidths_m = c(500, 1000, 2000, 4000),
    candidates = c(
      "geographic_only",
      "geographic_times_interpretable_environmental_similarity",
      "geographic_times_alphaearth_similarity"
    ),
    alphaearth = list(
      source_year = 2024,
      representation = "full_64_dimensions",
      standardization = "L2_normalize_nonzero_vectors",
      distance = "cosine_or_equivalent_chord_distance",
      pca = "fallback_only_if_distance_concentration_or_compute_diagnostics_fail",
      kernel = "exp(-0.5 * (chord_distance / bandwidth)^2)",
      primary_bandwidth_rule = "median_nearest_neighbour_distance_among_unique_external_locations",
      primary_bandwidth_observed = 0.232508650335214,
      bandwidth_status = "outcome_blind_external_covariate_distribution__no_carbon_response_used"
    ),
    interpretable_environment = list(
      variables = c(
        "NDVI_summer",
        "DW_floodedveg_prob",
        "DW_tree_prob",
        "relative_elevation_250m",
        "slope_deg"
      ),
      rationale = c(
        NDVI_summer = "summer vegetation productivity and cover",
        DW_floodedveg_prob = "surface wetland and inundated vegetation context",
        DW_tree_prob = "woody canopy context",
        relative_elevation_250m = "local landscape position",
        slope_deg = "local terrain gradient"
      ),
      scaling = "external_unique_location_median_and_MAD__no_carbon_response_used",
      distance = "per_variable_absolute_difference_on_external_robust_scales",
      kernel = "equal_mean_of_exp(-0.5 * standardized_variable_difference^2)",
      bandwidth_rule = "one_external_MAD_per_variable__no_multivariate_bandwidth",
      minimum_equivalent_core_weight = 0.1,
      variable_selection = paste(
        "prespecified mechanistic nonredundant set; omit spring NDMI and TCW",
        "because three community cores lack those seasonal values"
      ),
      amendment = paste(
        "The initial multivariate nearest-neighbour kernel was rejected by an",
        "outcome-blind scale diagnostic because flooded-vegetation probability",
        "alone reached 27 external MAD units and collapsed all weights. Equal",
        "averaging bounds each variable's contribution without using carbon."
      )
    ),
    precision = "retain_existing_conditional_proxy_without_ESS_or_Voronoi_adjustment",
    tuning_rule = paste(
      "Treat bandwidths as prespecified sensitivity cases; do not select a",
      "winner from eight local cores and claim independent validation."
    )
  ),

  hydrology_experiment = list(
    status = "supplemental_module__does_not_require_rerunning_the_119_band_export",
    result_status = "processing_passed__raw_strict_TWI_not_admitted",
    audited_on = "2026-09-25",
    admission_reason = paste(
      "conditioned elevation agreed but specific catchment area and TWI were",
      "not stable across DEMs, especially at the eight community cores"
    ),
    initial_dem_candidates = c("NRCan_CDEM", "Copernicus_GLO30_2024_1"),
    processing_grid_m = 30,
    aoi_buffer_m = 20000,
    breach_search_distance_cells = 20,
    flow_method = "D_infinity",
    accumulation_output = "specific_catchment_area_unlogged",
    twi_definition = "ln(specific_catchment_area_m / tan(slope_degrees))",
    strict_twi_rule = "retain_zero_slope_NoData_for_QA__do_not_silently_floor_slope",
    high_resolution_rule = paste(
      "Use HRDEM or another high-resolution DTM only where verified coverage",
      "intersects the AOI; never describe a resampled coarse DEM as new detail."
    ),
    processing = c(
      "buffer_the_AOI_for_flow_routing",
      "condition_the_DEM_and_resolve_flats",
      "use_a_named_D_infinity_or_multiple_flow_direction_method",
      "derive_specific_catchment_area_slope_and_true_TWI",
      "crop_to_an_interior_validity_area_and_retain_QA_layers"
    ),
    teaching_goal = paste(
      "Use WhiteboxTools through small explicit commands, keep intermediate",
      "rasters, explain each operation, and verify each result visually."
    )
  ),

  external_model_experiment = list(
    result_status = "completed__local_transfer_failed__not_admitted",
    completed_on = "2026-09-25",
    training_sources = c("RILEY", "CANPEAT"),
    first_target = "complete_0_30cm_carbon_stock_kg_m2",
    expected_training_profiles = 1512,
    final_local_transfer_test = "six_complete_community_0_30cm_profiles",
    models = c("glmnet_elastic_net", "ranger_random_forest"),
    optional_challengers = c(
      "small_regularized_DNN",
      "nonnegative_stack_from_identical_spatial_source_aware_out_of_fold_predictions"
    ),
    validation = "spatial_blocks_plus_source_study_grouping_where_available",
    community_role = paste(
      "untouched_final_transfer_test; local residual update only after a",
      "successful transfer gate (not reached in this run)"
    ),
    fitted_predictors = 105,
    external_unique_coordinates = 998,
    validation_groups = 180,
    external_random_forest_unique_location_mae_kg_m2 = 3.598431,
    external_mean_unique_location_mae_kg_m2 = 3.816128,
    community_external_mean_mae_kg_m2 = 4.529551,
    community_elastic_net_mae_kg_m2 = 4.611668,
    community_random_forest_mae_kg_m2 = 5.779218,
    community_random_forest_bias_kg_m2 = 5.276131,
    adoption_decision = paste(
      "reject external model as local prior: weak external signal and worse",
      "errors on the untouched six-profile community transfer test"
    ),
    prediction_extent = "community_AOI_only__never_a_Canada_wide_carbon_map",
    anisotropy = paste(
      "defer_until_hydrology_is_validated_and_residuals_show_directional",
      "structure; prefer DEM_defined_flow_or_cost_distance to estimating",
      "free anisotropy parameters from eight community cores"
    )
  ),

  adoption_gate = c(
    "improves_relevant_held_out_errors_without_hiding_per_core_failures",
    "stable_across_prespecified_bandwidth_sensitivity_cases",
    "defensible_landscape_interpretation",
    "acceptable_covariate_coverage_and_environmental_support",
    "no_unexplained_spatial_discontinuities",
    "baseline_retained_if_evidence_is_not_consistent"
  )
)
