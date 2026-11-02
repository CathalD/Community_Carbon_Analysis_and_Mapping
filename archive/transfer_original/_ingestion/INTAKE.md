# Transfer Oct 29 intake

Ingested 30 October 2026. This is a record of the supplied workflow and test data, prepared before the planning discussion. The original scripts and data are unchanged. The workflow has been read and its inputs inspected. Scripts 1–5 have not been executed as an end-to-end analysis.

## Source and local copy

- Repository: https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping
- Branch: `Transfer_Oct`
- Commit: `8533ac27a098e0e0d1e09b560cabe0038c80c7f3`
- Repository folder: `Transfer_OCT_29`
- Local folder: `/Users/cathaldoherty/Desktop/Community_Carbon_Analysis_and_Mapping-Transfer_Oct/Transfer_Oct29`
- Original folder contents: **30 files, 4,811,092 bytes**. The initial chat count of 31 included the folder entry itself.
- Contents: six R files, eleven TIFFs, four GeoPackages, four PNG maps, three CSV tables, one GeoJSON boundary, and a blank README.

`source_manifest.csv` records the original path, local path, byte count, Git blob SHA-1, and SHA-256 for every source file. Every original file was checked against the pinned Git commit. `source.json` records provenance for the transfer and supplemental files.

The two supplemental files under `reference_inputs/` come from elsewhere on the same branch:

1. `community_carbon_outputs_updated.xlsm`, from `Transfer_Oct/`, is the workbook named by Script 1. Its six simulated cores match the mapping CSVs.
2. `GEE_Raster_Exports.ipynb`, from the repository root, contains the covariate preparation and export workflow.

The user's attached overview is preserved verbatim as `Script_6_MRV_overview.txt`. Supplemental material and intake notes are separate from the 30 original files.

## Workflow as supplied

| Component | Inputs and role | Outputs / current state |
| --- | --- | --- |
| `Script_1.R` | Reads the Core Log and Peat Data workbook sheets, plots SOC and bulk density profiles, and sums section stocks by depth. Joins core coordinates. | Produces mapping CSVs with 0–30, 0–50, 0–100 cm and sampled-profile stocks. Uses `data/input/` and `Draft1/data/` paths. |
| `Script_2.R` | Reads mapping table and AOI, transforms to EPSG:3978, calculates area-mean extrapolations and 1 km IDW maps. Tunes IDW power and neighbour count by leave-one-core-out RMSE for each target. | Produces area polygons, depth-specific IDW rasters, tuning tables, maps, and clipped-area carbon totals. Uses both `Draft1/outputs/` and `outputs/`. |
| `Map_helper.R` | Defines context-map helpers using `data/context/north_america.gpkg`. | That context database is absent from this transfer. Script 2 sources the helper file, but does not call its context-map functions. |
| `Scrip_3.R` | Extracts covariates, checks coverage, fits small linear models, validates by leaving one core out, and maps predictions and extrapolation flags. | The filename really is `Scrip_3.R`. Its seven prescribed model formulas require predictor bands absent from the supplied elevation-only TIFF. With the supplied bands, those models would all be skipped. |
| `Script_4.R` | Compares field/area means, IDW, and a McMaster prior. Contains map differences and descriptive summaries. | Depends on earlier session objects and a Windows `Y:/` source path. Already-clipped prior TIFFs are supplied locally. Its IDW interval uses pixel SD divided by the square root of pixel count. That calculation is not a validated estimate of spatially aggregated prediction uncertainty. |
| `Script_5.R` | Starts with uncertainty comparisons, then defines a conditional log-Gaussian spatial update using a prior mean/SD and core observations. Includes held-out sampling-group predictions, predictive coverage, map exports and assumption records. | An unfinished draft with inconsistent object names and missing initialization. Its spatial update already assumes exponential covariance, a 5 km e-folding distance, and log-scale observation SD 0.15. These are fixed settings, not estimated parameters. |
| `GEE_Raster_Exports.ipynb` | Colab/Earth Engine preparation for terrain, Sentinel-2 vegetation, Dynamic World cover, Sentinel-1 radar, surface water, climate and observation-count QA. Includes export settings and band manifest generation. | Needs Earth Engine authentication and project access to run. Current code sets 30 m spacing and 2019–2023 dates; the supplied TIFF is a different, 250 m elevation-only product. The export loop lists all seven groups even though the earlier `covariate_groups` setting lists topography only. |

## Test data actually present

- The three CSVs (`mapping_data.csv`, `cores_for_mapping.csv`, `peat_cores_for_mapping.csv`) are byte-identical aliases of the same six-row, 30-column dataset. They do not constitute three independent datasets.
- Core and plot IDs are unique: `AOI-SIM-C01`–`C06` and `AOI-SIM-P01`–`P06`. Notes explicitly identify simulated measurements.
- The supplemental workbook contains 44 populated depth sections. Recalculating the Script 1 depth sums from its cached section-stock values reproduces all supplied CSV depth stocks to the exported two decimals. Cached workbook values were read; formulas were not recalculated or modified. Script 1's `skip=2` correctly finds the workbook headers with the available readxl installation.
- Maximum sampled depths are 30, 45, 60, 80, 100 and 120 cm. Six profiles reach 30 cm, four reach 50 cm, and two reach 100 cm. Current deeper stock columns retain shallower sampled totals, so their labels do not establish complete measurements to the target depth. None of the cores establishes mineral contact.
- Stocks are in **kg C/m²**. The overview's Mg C/ha examples require a factor of 10 conversion; total tonnes C are calculated from kg C/m² × area m² / 1,000.
- The AOI geometry is valid and contains all six cores. Its projected area is approximately **13,672.61 km²** in EPSG:3978.
- `topography.tif` has one band, `elevation_m`, at 250 m spacing, 573 columns × 536 rows. Elevation values are available at all six core locations. Relative elevation, slope, vegetation, cover, radar, water, climate, QA bands, `band_manifest.csv` and `settings.json` are absent from the supplied folder.
- The IDW and prior TIFFs align on a 1 km grid, 143 columns × 134 rows. Four `carbon_stock_*` TIFFs duplicate their corresponding `IDW_*` TIFFs byte for byte.
- `prior_MM_0_100.tif` and `prior_MM_sd.tif` have positive values at every core. Their stored band descriptions reference a 250 m McMaster source, but these supplied copies are already on the 1 km output grid. The original uncertainty definition and resampling provenance still need confirmation.
- `area_based.gpkg` contains the area polygons by depth, including both `area_based_0-100` and `area_based_0_100` layer names. `aoi_polygon.gpkg` stores the boundary. `IDW_based.gpkg` and `community_carbon_mapping.gpkg` contain raster coverage, rather than core point layers.
- The four PNGs were visually inspected. They show area-mean maps and core points at the four named depths, with a shared carbon scale.

## Execution readiness observations

All six R files pass R syntax parsing. This does not establish successful execution. Specific issues retained for planning:

1. The flat transferred files do not satisfy the scripts' expected `Draft1/data`, `Draft1/outputs`, `data/input`, and covariate-folder paths.
2. Script 4 and Script 5 depend on previous in-memory objects such as `area_based_gpkg` and `core_proj`/`core_sf`. Script 5 also references `prior_sd`, `rmse_100`, `prior_mean_file`, `prior_uncertainty_file`, `excluded` and `input_paths` without consistent initialization.
3. Script 5 line 107 has an unfinished assignment. R parses it together with line 109 as `values(idw_uncertainty) <- idw_cv <- (idw_uncertainty / idw_0_100) * 100`. It passes parsing but does not specify a usable IDW uncertainty model.
4. Script 5 mixes an `sf` AOI expectation with `SpatVector` AOI objects. It also converts field points but later extracts using a separate `core_sf` object.
5. Script 3 has no supplied bands matching its prescribed formulas. The available elevation band could support a separately specified model, but no model was changed during intake.
6. The shell's default R is 4.5.2 and lacks several required packages. A separate installation at `/Library/Frameworks/R.framework/Versions/4.6/Resources/bin/Rscript` is R 4.6.0 and has sf, terra, readr, dplyr, ggplot2, readxl, janitor, gstat, patchwork and jsonlite. Read-only checks used that installation. No packages were installed.

`raster_inventory.json` records raster dimensions, grid transforms, CRS definitions, bands, NoData and computed ranges. `read_only_audit.txt` records workbook ingestion, profile counts, geometry checks, core covariate/prior extraction, R parsing and software versions.

## Script 6 overview ingested

The attached document proposes Bayesian DSM at three depth targets, held-out validation, calibrated prediction intervals, pixel uncertainty maps, AOI totals and uncertainty, residual spatial analysis, sampling-sufficiency curves, adaptive sampling, and an MRV assessment report. Repeat-campaign change detection is proposed for a later Script 7. Script 6 code is not present in the source folder.

The overview's descriptions of the Wadoux papers and its performance targets are recorded as proposed methodological guidance. The referenced papers were not supplied or independently reviewed during this ingestion.

Points to retain for the forthcoming planning discussion:

- Match the target stock depth to actual observed support before comparing or fitting models.
- Distinguish posterior uncertainty in the latent stock from uncertainty in a future observed core and from uncertainty in the AOI total.
- Script 5 already includes assumed spatial covariance. Any residual covariance correction needs a defined role to avoid counting uncertainty twice.
- Six available cores can support only subset experiments up to six distinct cores. Larger sample-count curves require more observations or an explicitly defined simulation or forecasting design. There is no empirical evidence here yet for a 30-core threshold or a specified uncertainty reduction.
- The synthetic measurements permit workflow checks. They do not validate relationships between real terrain and measured soil carbon, establish MRV readiness, or demonstrate repeat-campaign change detection.

No workflow repairs, new models, new scientific estimates, or implementation plan were created during this intake.
