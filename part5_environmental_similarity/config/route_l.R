# Route L: variable organic-horizon carbon plus mineral carbon below the contact.
list(
  scenario = "route_l",
  analysis_title = "Variable-depth organic-horizon plus mineral soil carbon",
  mineral_depth_cm = 30,
  formal_precision_claim = FALSE,
  li_depth_cap_cm = 600,
  closure_tolerance_kg_m2 = 1e-4,
  li_source_class_definition = c(
    `1` = "Complete depth/carbon and positive spreads; curve eligible",
    `2` = "Usable depth/carbon; one or more spreads missing or nonpositive",
    `3` = "Usable depth/carbon below curve support (takes precedence over class 2)",
    `4` = "No usable depth/carbon pair"),
  community_change_definition = c(
    depth_effect = "Li curve at updated depth minus original Li carbon",
    curve_effect = "Updated organic carbon minus Li curve at updated depth",
    total_organic_change = "Updated organic carbon minus original Li carbon"),
  sensitivity = list(thickness_range_m = c(1000, 2000, 4000),
    observation_spread_cm = c(5, 10, 20), tail_transition_cm = c(25, 50, 100)),
  # The output uses Li's native 30 m grid as a carrier grid. It does not
  # claim that the cores establish 30 m-scale certainty.
  map_cellsize_m = 30,
  organic_layer_range_m = 2000,
  organic_layer_observation_spread_cm = 5,
  # First-pass spatial scale for organic carbon-density profiles. Kept separate
  # from thickness and exposed for future local/regional covariate refinement.
  organic_carbon_profile_range_m = 2000,
  # Provisional first-pass mineral mapping scale. This is deliberately kept
  # separate from the organic-thickness setting and must be sensitivity-tested.
  mineral_carbon_range_m = 2000,
  organic_layer_label = "measured organic-layer thickness",
  organic_contact_rule = "field_observed_preferred__core_layer_classification_fallback",
  organic_soc_threshold_g_kg = 80,
  paths = list(
    aoi = "data/raw/aoi.geojson",
    core_layers = "data/raw/core_layers.csv",
    core_locations = "data/raw/core_locations.csv",
    field_contacts = "data/raw/organic_mineral_contacts.csv",
    li_organic_thickness = "data/external/McMaster_PeatDepth_30m_cm_version1.0.tif",
    li_organic_thickness_uncertainty = "data/external/McMaster_Uncertainty_PeatDepth_30m_cm_version1.0.tif",
    li_organic_carbon = "data/external/McMaster_Carbon_30m_kgm2_version1.0.tif",
    li_organic_carbon_uncertainty = "data/external/McMaster_Uncertainty_Carbon_30m_kgm2_version1.0.tif"
  ),
  li_uncertainty_definition = "ensemble_spread__not_calibrated_predictive_sd",
  li_weighting = "inverse_spread_squared__reported_as_sensitivity",
  li_carbon_units = "kg_C_m2__whole_peat_profile",
  li_carbon_depth_relationship_published = "C_kg_m2 = 0.38 * peat_depth_cm + 17.14",
  # Estimated from the supplied co-registered rasters within the AOI. Retain
  # separately from the paper's rounded published equation above.
  li_carbon_depth_slope_raster_matched = 0.375097007,
  li_carbon_depth_intercept_raster_matched = 17.139840974,
  li_carbon_uncertainty_definition = "plus_minus_1_SD__base_model_disagreement__not_independent_of_depth_uncertainty",
  li_carbon_role = "paired_depth_derived_prior_and_consistency_check__not_a_second_independent_prior",
  li_carbon_curve_rule = "linear_cumulative_curve_through_zero_using_Li_mean_density__only_where_Li_depth_is_at_least_40_cm",
  li_carbon_minimum_depth_for_curve_cm = 40,
  organic_carbon_tail_transition_cm = 50,
  community_curve_tail_rule = "preserve_observed_core_curve__then_exponentially_decay_cumulative_residual_to_local_Li_curve",
  community_carbon_update_rule = "Li_baseline_at_blended_depth_plus_thickness_core_influence_times_distance_weighted_core_residual",
  notes = paste(
    "No peat-presence model. Organic thickness is Li plus observed contacts;",
    "mineral carbon is spatially extrapolated from complete measured mineral layers only;",
    "Li carbon is depth-derived, so it is not fused independently from Li depth;",
    "organic carbon preserves measured core curves and uses a Li-informed tail beyond measured depth;",
    "mineral carbon is provisional; support/extrapolation diagnostics are not uncertainty estimates; Li-incomplete cells are explicitly flagged."
  )
)
