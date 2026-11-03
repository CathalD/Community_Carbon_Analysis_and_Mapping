# Fort Severn — full workshop comparison

> **GitHub snapshot:** large local inputs and generated spatial/model files are omitted. Current figures and compact results are included. See [package contents and reproduction instructions](docs/GITHUB_PACKAGE.md) and the [file manifest](docs/github_package_manifest.csv).

Open `FS_test.Rproj` in RStudio. This project applies the workshop to Fort Severn with **four method branches and four corresponding prior updates**: raw cores, area average, IDW and covariate regression. Sothe has four stock updates. Li has four depth updates and four organic-carbon updates.

## Run

```r
source("run_workflow.R")
run_fs_workflow()
```

To rerun all model stages using the already-prepared inputs:

```r
run_fs_workflow(c("1", "2", "3", "4", "5", "5c", "3d", "5d", "verify", "results"))
```

In this computer's terminal use `/Library/Frameworks/R.framework/Resources/bin/Rscript run_workflow.R 1 2 3 4 5 5c 3d 5d verify results`. The Homebrew R installation lacks spatial packages. Required packages: terra, sf, readr, dplyr, ggplot2, jsonlite, gstat and patchwork. Logs and status are under `outputs/`.

| Step | Script / purpose |
|---|---|
| prepare | `prepare_inputs.R`: crop original priors, prepare 119 environmental bands and native core extractions. |
| 1 | `script_1.R`: workshop field calculations and reference comparisons, all eight FS cores. |
| 2 | `script_2.R`: unchanged workshop area averages and IDW tuning/maps/totals. |
| 3 | `script_3.R`: covariate regressions, holdouts and AICc selection for the Sothe reporting group. |
| 4 | `script_4.R`: align and compare local maps and Sothe prior. |
| 5 | `script_5.R`: Sothe raw-core spatial, area-based and IDW updates. |
| 5c | `script_5c_sothe_regression_update.R`: fourth Sothe update using the selected regression map. |
| 3d | `script_3d_li_model_suite.R`: four local methods for Li-matched depth and carbon, plus cumulative-carbon curves. |
| 5d | `script_5d_li_update_suite.R`: four Li updates per variable, nested core holdouts and sensitivity. |
| verify | `verify_outputs.R` / `verify_full_suite.R`: completeness, data lineage, zeros and independent numerical checks. |
| results | `show_full_suite_results.R`: current comparisons of prior plus four updates and all four local methods. |

## All eight cores and the workshop depth rule

The workshop's `stock_to_depth()` function is unchanged: integrate measured sections to **min(target depth, sampled bottom)**, prorating a section crossing the target. Shallow cores remain included. No carbon below the sampled bottom is invented. Depth support is in `outputs/script1/depth_coverage.csv`. SOC is converted from g/kg to workshop percent; stock is independently checked from raw density, carbon and thickness. Modeling retains full precision.

The workshop reference data keep their original eligibility rules. Their depth support can differ from these field groups. The deeper 50/100 cm field columns remain measured-to-bottom products, not complete fixed-depth inventories.

## Sothe: four updates

Outputs are under `outputs/script5/sothe/stock_0_30/`:

- `raw_cores`: workshop spatial log-Gaussian update.
- `area_based`: pointwise fusion with the area average and between-core SD.
- `idw`: pointwise fusion with tuned IDW and working held-out prediction error.
- `regression`: pointwise fusion with the AICc-selected Script 3 map and fixed-model holdout error.

The mean prior is Sothe S30. S100 full Q95−Q05 width is converted with `width/(2*qnorm(.95))` and used directly as the requested absolute S30 SD proxy, without scaling by 0.3. This is a cross-depth sensitivity, not calibrated S30 uncertainty. All four maps use the same 3,873 prior cells. They are alternatives using the same observations; do not combine them as independent evidence.

FS-04's containing Sothe pixel is NoData. Its conditioning prior uses the nearest valid paired source pixel, 129.72 m away, within a 500 m limit. The source map mask is unchanged and extraction is audited. Script 3 selects Mean_only, so the fourth branch legitimately retains an intercept-only regression.

## Li: four methods and four updates per variable

The user selected **all eight organic-contact estimates**, retaining evidence labels. Two contacts are field-confirmed. Six are reproduced from the current project's layer rule: top of the first layer below 80 g SOC/kg, with confirmed contacts taking precedence. No Route L mapping, blends or tails are called.

| Core(s) | Organic depth, cm | Evidence |
|---|---:|---|
| PM-01-A | 14.50 | Confirmed |
| PM-03-A | 30.73 | Confirmed |
| PM-02-A | 15.20 | Layer-classified |
| FS-05 | 10.00 | Layer-classified |
| FS-01, FS-02, FS-03, FS-04 | 0 | Layer-classified zero organic thickness |

Layer classifications are assumptions about contacts, not independent confirmation of peat presence/absence. Labels are retained in every core table. Organic stocks integrate measured layers to those boundaries. Zero-depth responses have zero organic stock; their mineral carbon remains separate from Li's organic pool.

Local methods are under `outputs/script3d_li_suite/{depth,carbon}/`: raw observations/points, area mean, tuned IDW and AICc-selected regression. Raw observations remain points, not an invented interpolated surface. All 119 covariates are audited; a small prespecified candidate set is fitted to avoid 119 coefficients for eight cores.

Updated maps are under `outputs/script5/li_suite/{depth,carbon}/{raw_cores,area_based,idw,regression}/`. `completion_manifest.csv` enumerates all eight updates; `suite_summary.csv` records means, validation and support counts.

### Explicit scale adaptation for zeros

The workshop lognormal model cannot take log(0). To retain all authorized labels, **Li applies the workshop Gaussian conditioning/precision algebra on the original measurement scale**. This is an explicitly labelled Gaussian sensitivity, not an unchanged lognormal carbon model. No offset is added and no zero response is dropped.

Area/IDW/regression fusion uses between-core SD or held-out model RMSE as working errors. Spatial depth uses an assumed 5 cm observation spread with 2.5/5/10 cm sensitivity; carbon retains zero measurement error. Spatial range sensitivity is 1/2/4 km.

Negative local regression predictions are flagged and not used as physical evidence: the prior is retained there with `update_applied=0`. Negative Gaussian posterior means remain visible with flags and `probability_below_zero`; they are not silently clipped or interpreted as physical negative peat/carbon. The current run has 33 such depth-GP cells and 14 carbon-GP cells.

FS-02/FS-04 have native prior NoData; paired-source conditioning priors are taken 21.85/18.96 m away and audited. Valid zero-valued depth priors remain valid. Source rasters and masks are unchanged.

Outer holdouts refit local methods, retune IDW, reselect regression and estimate errors using only training cores. These remain eight nearby, partly classified observations. Li ensemble spreads are not calibrated prediction errors. Depth and carbon remain dependent marginal comparisons, not independent posteriors to combine into a joint inventory.

### Carbon curves and depth-to-carbon diagnostics

`outputs/script3d_li_suite/carbon_curve/` contains measured cumulative curves and a through-zero fitted carbon-depth curve. Only the four positive organic profiles contribute slope information/residual degrees of freedom. Zero profiles are shown at the origin. Maximum measured support is 30.73 cm.

Each updated depth branch includes `curve_derived_carbon.tif`: the local curve evaluated at updated mean depth only within measured support, plus a support flag. This connects each depth update to the curve without inventing a deep tail. It is distinct from direct organic-stock updates and does not propagate joint depth/curve uncertainty. Temporal accumulation rates require chronology.

## Current results and additional data

Current figures are in `outputs/results/full_suite/`; see `docs/run_results.md`. More confirmed contacts and deeper complete peat profiles would materially improve Li. `data/field/additional_contact_template.csv` lists six contacts needing confirmation. Field/lab error estimates would improve the working likelihoods.

`archive/workshop_source/` preserves original workshop scripts and hashes. Earlier two-contact scripts/results (`3b`, `5b`, `outputs/script5/li/`) remain historical and are excluded from the current default run. Earlier documentation is in `archive/previous_two_contact_run/`. No original workshop or production Route L analysis is modified. Set source paths in `config.R` when rebuilding inputs after moving the folder.
