# Part 5 — Environmental similarity and external donor observations

This self-contained folder extracts the relevant data, code, provenance and prior findings from **CommunityCarbonMap_v3** for section **5 of the supplied methods description**. It is separate from the workshop's existing `script_5.R` (Bayesian map updates). Open the workshop RStudio project as usual; the commands below locate this folder and restore your working directory afterward.

## Run the donor comparison

From the workshop RStudio console:

```r
source("part5_environmental_similarity/run_part5.R")
```

Or open a terminal in this folder and run `Rscript run_part5.R`.

Only base R is needed for this starter run. It uses the copied tables; no account, download, installation, original-project path, or raster processing is needed. A rerun replaces this folder's generated comparisons.

The run first reproduces the original arithmetic-stock donor comparison, then runs a new, explicitly exploratory log-carbon weighted-donor comparison following the supplied section 5 equation. It compares source log mean, a Gaussian environmental kernel, and an AlphaEarth kernel separately for CanPEAT and Riley. It reports external grouped holdouts, whole-database transfer, and the original community transfer test. Settings are in `part5_settings.R`.

Inspect:

- `outputs/part5_log_donors/metrics.csv`: MAE, RMSE, bias and log RMSE by method/source/test.
- `outputs/part5_log_donors/predictions.csv`: every prediction, environmental support, donor spread and effective donor count.
- `outputs/part5_log_donors/community_top_donors.csv`: donor identities, source studies, ecosystems, coordinates and normalized weights.
- `outputs/part5_log_donors/training_pool_scaling.csv`: donor-only scales and AlphaEarth bandwidth for each test pool.
- `outputs/part5_log_donors/community_comparison.png`: original-study measured stocks versus geometric stock estimates.
- `outputs/covariate_analysis/part15_external_similarity_transfer/`: rerun of the unchanged historical 25/50-neighbour raw-stock comparison.

The log run estimates `sum(w_i * log(C_i))` with normalized similarities. Its exponent is a **geometric stock estimate**, not a lognormal arithmetic predictive mean. No fitted predictive distribution, GP covariance, interval coverage or posterior precision is claimed. Donor log spread is descriptive variation among donor outcomes; effective donor count describes concentrated weights. Neither is a prediction-error variance.

## Included observations

| Table | Contents |
| --- | --- |
| `data/derived/covariate_analysis/external_core_covariates_riley_canpeat.csv` | 1,564 unique donor records: 1,217 CanPEAT and 347 Riley, carbon responses, provenance, 119 predictors and QA |
| `data/derived/covariate_analysis/external_training_0_30cm.csv` | Existing 1,512 complete 0–30 cm training records |
| `data/derived/covariate_analysis/community_core_covariates.csv` | Eight **original-study** cores, six with complete measured 0–30 cm stocks |
| `data/derived/external_carbon_profiles/` | Harmonized Riley/CanPEAT cores and layers, source manifest, wider stock table and Canadian point-export inventory |
| `data/external/Open_Core_Data/` | Original local Riley workbook/CSV and CanPEAT core/profile extracts |
| `data/external/Hierarchical_Data/` | Twelve newer multiscale point shards: 1,564 donors plus eight original community targets |
| `data/raw/gee_imports/` | Historical point exports used to prepare the 119-predictor tables |

The common-covariate comparison retains **1,504 profiles at 996 unique source/coordinate sites**, consisting of 855 CanPEAT and 141 Riley sites. Missing predictors explain the difference from 1,512 eligible responses. Profiles sharing a source/coordinate are averaged to one arithmetic stock before logging or neighbour selection; their first exported covariate vector represents that site, as in the historical comparison. Do not interpret repeated profiles as independent locations. The wider stock database includes other collections for context; the executable donor comparisons use Riley and CanPEAT only.

`docs/data_inventory.csv` describes table sizes and roles. `docs/prepared_table_columns.csv` inventories the prepared columns. `source_manifest.csv` records original absolute paths, copied-file hashes and sizes. `bundle_checksums.csv` records the delivered files.

## Predictors and validation

The interpretable environmental set is NDVI_summer, DW_floodedveg_prob, DW_tree_prob, relative_elevation_250m and slope_deg. AlphaEarth uses all 64 dimensions, L2-normalized. The existing export also retains seasonal vegetation/moisture indices, terrain, neighbourhood summaries and observation-quality counts. QA counts remain separate from predictors. Climate and substrate are potential future families, not claims about available bands in this first export. A later 47-band multiscale recipe and its point extracts are included separately.

In the new log run, environmental medians/MADs (SD fallback) and the AlphaEarth nearest-neighbour bandwidth are fitted from each training pool only. The environmental kernel is `exp(-0.5 * sum((z_target-z_donor)^2 / ell^2))`, with fixed `ell = sqrt(5)` for five variables. These are exploratory settings, not selected as optimal. No local response selects settings.

External holdouts reuse the original five connected groups built from 100 km blocks, source studies and rounded coordinates. Shared groups cannot cross folds. The split was originally response-balanced; preserve it as the historical evaluation design. Whole-database transfer is also reported. The original community outcomes have been examined in earlier experiments, so they are a held-out transfer example for this run, **not a fresh independent test for future model development**.

## Using the workshop's simulated cores or a new AOI

The workshop's six synthetic cores are a different dataset from the eight original-study cores supplied here. This folder does **not** substitute old covariates for the synthetic locations or join their outcomes by row order. Script 1's carbon calculations and capped-at-recorded-maximum depth convention remain unchanged.

A comparable full 0–30 cm external response and a synthetic core capped above 30 cm are different targets. Keep the workshop convention explicit and review target compatibility before evaluating or conditioning a donor model. The included original-study transfer example uses the original study's complete-depth eligibility, not a new workshop exclusion rule.

For new local cores and prediction cells, extract exactly the donor predictors with the same source versions, years, seasonal definitions, QA, neighbourhood sizes and spatial support. The copied `scripts/gee/Covariate_HBJB.ipynb`, `route_l_multiscale_covariates.py`, `canada_external_core_sampling.py` and `aoi_community_core_sampling.py` preserve the first recipe; the donor sampler and micro-pilot files preserve the later multiscale recipe. They contain historical cloud project/asset settings that must be replaced with accessible ones. They are reference extraction code; no remote exports have been started. Do not mix the two predictor recipes or treat the workshop's topography/water/optical files as equivalent just because names resemble each other.

After matching the new observations and map cells to the donor recipe, build a separate target table keyed by stable IDs. The present starter deliberately reproduces the original targets; editing only `community_core_covariates.csv` is insufficient because the historical script asserts its six-core study design. A new-target adapter or model is a subsequent implementation, and must not bypass these compatibility checks.

## Optional historical multi-depth example

```r
source("part5_environmental_similarity/run_hierarchical_example.R")
```

This runs the copied multiscale point pipeline, its own self-test, measured shallow validation and source-mean baselines. It uses five intervals (0–15, 15–30, 30–50, 50–100, 100–200 cm), local/regional/broad covariate blocks and depth exponents. Results go to `outputs/covariate_analysis/part16_hierarchical_similarity_micro_pilot/fs_aoi_point_test_v1/`. It is a raw-stock historical comparison, not the log model above. Raster mapping is optional in the underlying script and requires `terra` plus a matching new raster; no AOI map is generated by this wrapper.

## Existing findings and further modelling

`reference_results/` preserves the external elastic-net/random-forest and donor-comparison results separately from regenerated outputs. The random forest had external signal but transferred poorly to the original community cores; similarity methods showed source bias and weak local support. The later hierarchical point example likewise did not justify mapping the larger AOI. See `docs/source_notes/external_model_part5_findings.md`, the preserved Part 15 comparison, and `ANALYSIS_REPORT.md`. These are prior study findings, not measured performance on the workshop's synthetic data. The original `reference_results/part15_external_similarity_transfer/EXTERNAL_COMPARISON.md` contains image links to the accompanying copied figures.

The copied `scripts/56_prepare_external_carbon_model.R`, `57_fit_external_carbon_models.R` and `R/covariates/external_carbon_model.R` retain the original model experiment. They require `sf`, `Matrix`, `glmnet` and `ranger`. The original ingestion/preparation scripts are retained as provenance and may expect additional original-project inputs; use the included prepared tables for the self-contained run. No packages are installed automatically.

A GP fit, joint local/map covariance, conditioning on new local measurements, lognormal predictive moments, and map exports belong to sections 4/7 of the supplied description and are **not implemented or implied by this data extraction**. The donor database and validated examples here provide the section 5 starting point. A full GP requires a separate model specification and validation; do not assign donor spread as its likelihood variance.

## Source context

`docs/request_context.txt` preserves the supplied methods description. Source documentation is historical context, not instructions to modify the workshop. `config/route_l.R` is copied solely because the historical script checks its presence; this runner neither sources that file nor updates Route L. Original source attribution, measured/estimated flags and unresolved licence notes remain in the copied records. No source-data licence is changed by this local transfer.
