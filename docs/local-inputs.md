# Real local inputs and optional stages

Keep the complete example: `source("run_all.R")`. The shorter `run_my_data.R` still runs field calculations and IDW only. The additional `run_local_workshop.R` runs those calculations and reference comparisons, then runs optional ecosystem, regression and Bayesian stages when their input paths are supplied. Edit `local_settings` at the top of that file. No scientific inputs are created from an AOI alone.

## Begin with measurements and an AOI

Put measurements in `data/input/my_field_data.csv` with the [eight field columns](../data/input/field_columns.md); put the polygon in `data/aoi/my_aoi.geojson`. All local interfaces use soil organic carbon at 0–30 cm below the ground surface, in kg C/m². Missing surface or interval coverage gives no fixed-depth stock. Entered OC is authoritative. Bulk density and OC must refer to the same dry material and sample volume; resolve coarse-fragment corrections before import. West longitudes are negative. Projected AOIs need a correct embedded CRS; core coordinates are WGS84 degrees.

Run `source("run_local_workshop.R")` with later paths left `NULL`. It writes field results and reference comparisons under `outputs/my_data/part1/`, maps under `outputs/my_data/maps/`, and `local_grid.gpkg` plus `local_grid_template.csv` in `outputs/my_data/`. Review `local_stage_status.csv` for completed and unconfigured stages. This is a successful partial local analysis, not a full inventory.

The template has `cell_id`, `x_m`, `y_m`, `area_m2`. Save copies under `data/input/`, add the columns below, and enable their paths in `local_settings`. Retain all template rows and geometry columns unchanged. The code matches IDs and checks coordinates and clipped area, preventing accidental reuse after a grid/AOI change. Keep the resolution fixed between input preparation and rerun.

Raster pixel centres (`x_m`, `y_m`) are prediction locations. Boundary pixel centres may fall outside the clipped AOI. Values are point predictions displayed on clipped cells, not calibrated cell means; multiplying by clipped area is a numerical inventory approximation. GeoTIFF boundary pixels extend beyond the AOI; use `local_grid.gpkg` for the precise mask/area weights. Never infer exact inventory area by counting nonmissing raster pixels.

## Ecosystem totals

Set `ecosystem_cells` to `data/input/my_ecosystem_cells.csv`: the grid columns plus `ecosystem`. Obtain labels from a real land-cover/ecosystem product using an explicitly chosen rasterisation rule. Do not assign labels from nearest cores. Labels must match `ecosystem` in the field table. Every mapped ecosystem needs usable stock observations. At least one allows calculation, but does not establish representative sampling or quantify uncertainty.

The output `local_area_totals.csv` contains ecosystem areas, sample counts, means and tonnes C. `local_group_map.png` displays the assigned means. Summing the tonnes column gives the whole-AOI total. Mixed cells, unsampled classes, unequal sampling effort and inappropriate coarse-fragment conventions require scientific decisions outside this simple estimator. Choose a resolution suitable for your ecosystem map; this route uses one class per clipped cell.

## Landscape regression and held-out assessment

Set both `predictor_cells` and `model_observations`. The cell table adds numeric `wetness` and `elevation_m` to the grid columns. Supply predictors extracted from real maps at the prediction locations in the grid CRS; elevation is metres, and wetness must use the same defined index and scale at samples and prediction locations. Reproject and align source rasters deliberately in your GIS; these CSV interfaces do not silently resample rasters. Retain source, date, resolution, interpolation rule and units with your inputs.

The observation table requires `core_id,longitude,latitude,stock_0_30_kg_m2,wetness,elevation_m,split`. Stocks must be complete, comparable 0–30 cm observations; distinct locations lie inside the AOI. Assign `Training`, `Spatial holdout` or `Buffer` without inspecting held-out responses. Set a positive `holdout_buffer_m` appropriate to spatial dependence and the distance over which the model will predict. The code rejects points closer across the training/test divide than that distance. Minimum counts of ten training and five test locations are software safeguards, not proof of adequate sampling.

Outputs include `local_holdout.csv`, `local_model_metrics.csv`, `local_model_coefficients.csv` and `local_model_map.png`. The code assesses the training fit, then refits all supplied observations for mapping. Predictor selection or tuning using test responses would invalidate that assessment. Additional spatial folds and residual checks are needed for a real scientific claim. The final map is not used automatically as Bayesian evidence; the local update uses the IDW branch.

## Prior comparison and Bayesian update

Set `prior_cells` to the grid table with `prior_mean_kg_m2`, `prior_sd_kg_m2`, `new_sd_kg_m2`. Values must represent the same carbon pool, 0–30 cm surface depth and point support. Reconcile original map cell support before using this interface. An SD must describe error of the corresponding estimate, not just variation among cores or a regression prediction of a different quantity.

Review shared data, systematic bias and correlated errors. Set `independent_errors_confirmed=TRUE` only when the independent Gaussian error approximation is defensible. `local_bayesian.csv` retains prior, new and posterior means/variance/SD. `local_bayesian_map.png` displays the update. Compare `new_mean_kg_m2` with `prior_mean_kg_m2` before interpreting it. Negative Gaussian bounds are possible; no total uncertainty is inferred without cross-cell covariance.

## Ordinary kriging for real data

The maintained `scripts/04b_kriging.R` teaches covariance fitting on simulated observations. It deliberately does not auto-fit a four-core field variogram. Use `R/kriging.R` as a worked method, then design and review the real variogram, trend, support, anisotropy and spatial assessment before a field kriging analysis. The constrained exponential range search is a teaching choice, not an automatic model-selection prescription for every Canadian AOI.

## Reruns and limits

The full course scripts invalidate their named outputs and dependent stages before running. On a stage error they clear that stage’s partial results too. Learner notes are preserved. `outputs/run_complete.txt` exists only after `run_all.R` completes; opening an old image never recomputes it. Local runs clear their owned advanced products and earlier comparisons before work starts; they never fall back to demonstration inputs.

Automatic resolution targets roughly 2,500 bounding-box cells. A hard cap of 10,000 cells and two million distance pairs prevents unbounded allocation. Thin, irregular and multipart polygons use clipped areas. The AOI-centred equal-area projection works throughout Canada but distorts long-distance interpolation on continental extents; the code warns above a 1,000 km span. Use regional domains for scientific mapping. Do not override a geometry/CRS warning without checking coordinates.
