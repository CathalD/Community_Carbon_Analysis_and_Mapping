# Fort Severn — Carbon Analysis and Mapping workshop

Open `FS_test.Rproj` in RStudio. This project runs the supplied workshop with the eight measured Fort Severn cores from CommunityCarbonMap_v3. It does not run or source Route L models, blends, tails, or existing analysis outputs.

**Completion status:** the current Li implementation is partial. It includes local depth/carbon-curve fits and the raw-core spatial update applied separately to depth and carbon. It does not include a four-model/four-update Li comparison: area-average, IDW and covariate-regression Li branches and corresponding updates remain unimplemented. Passing the existing verification checks confirms the implemented calculations, not completion of these missing branches. The Sothe comparison contains the prior plus three update methods.

## Run

```r
source("run_workflow.R")
run_fs_workflow()
```

Or run individual steps, in order:

```r
run_fs_workflow("prepare")
run_fs_workflow("1")
run_fs_workflow("2")
run_fs_workflow("3")
run_fs_workflow("4")
run_fs_workflow("5")
run_fs_workflow("3b")
run_fs_workflow("5b")
run_fs_workflow("verify")
```

The runner uses a fresh R process for each script and saves logs in `outputs/logs`. In this computer's terminal use `/Library/Frameworks/R.framework/Resources/bin/Rscript run_workflow.R`; the Homebrew R installation does not have the spatial packages. Required packages: terra, sf, readr, dplyr, ggplot2, jsonlite, gstat, patchwork.

| Script | FS adaptation / outputs |
|---|---|
| `prepare_inputs.R` | Crop the original Sothe and Li rasters; prepare all 119 terrain, surface-context and AlphaEarth predictors on the workshop 1 km grid; native covariates at cores; input manifests. |
| `script_1.R` | Workshop profile/stock calculations and reference comparisons, using raw FS CSVs through a column-name/unit adapter. |
| `script_2.R` | Original workshop area averages, IDW tuning, maps and totals, unchanged. |
| `script_3.R` | Workshop linear models, whole-core holdouts, AICc selection and maps, using the expanded FS covariates. |
| `script_4.R` | Compare area, IDW, selected covariate model and Sothe S30 prior. |
| `script_5.R` | Workshop raw-core spatial Bayesian, area-based and IDW updates using Sothe S30 and the S100 uncertainty proxy. |
| `script_3b_peat_depth_and_curves.R` | Fit local peat-depth and cumulative-carbon models before the Li update. Two confirmed contacts support exploratory models only. |
| `script_5b_li_update.R` | Apply the workshop exponential log-Gaussian update separately to Li depth and carbon priors; holdouts and sensitivity. |
| `verify_outputs.R` | Check measured integration, all-eight-core inclusion, numerical conditioning, prior/posterior checks and agreement with the original workshop equations. |

## Depth calculations — all eight cores included

The active workshop's `stock_to_depth()` function is retained unchanged: integrate measured sections to **min(target depth, sampled bottom)**, prorating a section crossing the target. A core ending before 30 cm is included, with no invented carbon below its sampled bottom. This rule applies to the field 30, 50 and 100 cm columns. Depth coverage is saved in `outputs/script1/depth_coverage.csv`.

The workshop reference-data preparation uses additional depth eligibility rules. Its existing reference distributions are preserved, so their depth support can differ from the included shallow field profiles. The workshop README and reference rules do not override the requested field inclusion rule.

SOC is converted from raw g/kg to workshop percent, and stock is independently calculated as `SOC_gkg * BD_gcm3 * (1-coarse_fraction) * thickness_cm / 100` kg C/m². Full precision is retained for modeling rather than rounding to two decimals.

## Analysis 1: Sothe 0–30 cm

The mean prior is Sothe's 0–30 cm product. As requested, the **absolute** 0–100 cm full 90% interval width is used to estimate a working SD with `width/(2*qnorm(.95))`, matching workshop Script 5's conversion. It is not scaled by 0.3. The resulting S30 uncertainty, intervals and variance reduction are **conditional sensitivity estimates**, not calibrated S30 prediction intervals. Holdouts and half/double-SD and 1/2/4 km range sensitivity are saved alongside the raw-core update.

FS-04's containing Sothe pixel has NoData in both mean and width. To include this measured core, its conditioning prior uses the nearest valid paired source pixel, 129.72 m away, within a maximum allowed 500 m. The original raster's NoData mask is preserved. `core_prior_extraction_audit.csv` records the original and used values for all eight cores; this explicit FS input adaptation is the only addition to Script 5's prior extraction.

`outputs/script5/sothe/stock_0_30/` contains `raw_cores`, `area_based`, and `idw`. These are alternative workshop updates of the same measurements; do not combine them as independent evidence.

## Analysis 2: Li variable-depth peat carbon

All four native Li products are supplied: peat depth and spread (cm), whole-peat-profile organic carbon and spread (kg C/m²). The added local models use **explicit field-confirmed organic/mineral contacts**, currently PM-01-A (14.5 cm) and PM-03-A (30.73 cm). Sampled bottoms of other cores are not silently converted to peat depths. All eight remain in Analysis 1.

The depth baseline is the workshop's mean-only model. The cumulative carbon-depth curve is a pooled regression through zero, using one complete organic-stock endpoint per core. Measured section curves are also exported. Two contacts do not support a credible covariate peat-depth model or validated deep carbon curve. Mapping the local curve is restricted to measured depth support, so unsupported cells are explicitly NoData.

Li depth and carbon are updated **separately as marginal comparisons**, using the workshop spatial covariance equations. Carbon is depth-derived in Li's product; these results are dependent and must not be fused into a joint independent inventory. Source ensemble spreads are not calibrated predictive SDs. The depth run uses an explicit 5 cm observation-spread assumption with 2.5/5/10 cm sensitivity. Carbon retains the workshop zero-measurement-error baseline. Outputs are exploratory.

## Data to add for a stronger Li analysis

- Field-confirmed organic/mineral contact depths for the remaining six cores, or explicit peat absence. `data/field/additional_contact_template.csv` lists the missing core IDs; transfer confirmed entries into `data/field/organic_mineral_contacts.csv`. Confirmed zero peat needs a presence/absence treatment before the positive log-Gaussian update; it must not be changed to an arbitrary positive depth.
- More complete peat profiles, especially deeper than 30.73 cm, with measured bulk density and carbon. These support the variable-depth carbon curve and its extrapolation uncertainty.
- Field/laboratory error estimates if available. Current observation errors are working assumptions.
- If “carbon accumulation” means a rate over time, add chronology/age measurements; the current output is cumulative stock versus depth.

No additional inputs are required to rerun the current Sothe analysis or exploratory two-contact Li analysis.

## Covariates and provenance

All 119 downloaded environmental bands are retained and audited. The fit compares a small, prespecified set of one/two-predictor models rather than 119 coefficients fitted to eight cores. Predictors missing at any core are excluded without looking at the carbon response, retaining all eight cores for a fair comparison. AICc selects the model, as in the workshop; reported holdout scores are not nested estimates of model-selection performance.

Raster covariates are bilinearly sampled at centres of the 1 km analysis grid, not described as 1 km area averages. Field fitting uses native-raster extraction at core locations. Incomplete tiles stay missing, and prediction flags identify values outside observed predictor ranges. `data/covariates/band_manifest.csv` records all source tiles. Large inputs remain in their original locations; small cropped priors and predictor products are local. `prepare_inputs.R` needs the original source folders if rebuilt after moving this project; set paths in `config.R`.

`archive/workshop_source/` preserves the exact supplied scripts and README. `docs/workshop_source_manifest.csv` records their hashes. `docs/methods_and_changes.md` explains adaptations and source definitions. Existing projects and their outputs are untouched.
