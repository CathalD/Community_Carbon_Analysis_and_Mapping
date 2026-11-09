# CommunityCarbonMap v3 — Route L

Route L compares original Li/McMaster organic depth and carbon with community cores, then produces a locally updated organic-carbon map. A separate provisional mineral component covers 0–30 cm below the local organic/mineral contact. Changes are revisions to estimates, not measured carbon loss.

## Start in RStudio

Open `CommunityCarbonMap_v3.Rproj`, then run:

```r
source("run_route_l_step_by_step.R")
run_route_l_step("preflight")
run_route_l_step("test")
```

Loading the runner does not start the analysis. Run each subsequent command separately and inspect its outputs. The runner resolves the project from its source location or current directory. For an explicit location use `options(route_l.project_dir = "/path/to/CommunityCarbonMap_v3")` before sourcing it.

Required analysis packages are `terra` and `sf`. `renv.lock` records the working package versions without automatically activating a new library; see `docs/run_in_RStudio.md` before restoring. Raw data and external rasters must already be present. `config/route_l.R` defines assumptions and input paths. Exact core coordinates remain internal review material.

## Steps and review points

Use `run_route_l_step("ID")`. A step overwrites its own generated outputs. After changing inputs or assumptions, rerun dependent steps in order; the workflow does not automatically rebuild stale products.

| ID | Purpose and prerequisites | Inspect |
|---|---|---|
| 01 | Validate raw inputs and configuration | Run contract |
| 02 | Ingest/QA, after 01 | Core and layer audit |
| 03 | Harmonize depths/contacts, after 02 | Contacts and layer segments |
| 20 | Blend Li and measured organic depth, after 03 | Depth map, source flags, geographic influence |
| 21, 22 | Prepare/map provisional mineral evidence, after 03 | Completeness and extrapolation diagnostics |
| 23, 24, 25 | Review measured organic support and partial accounted stock, after 20–22 | Observed versus extrapolated support; not final total |
| 26a | Audit paired Li depth/carbon, after 20 | Li curve eligibility and product relationships |
| curves | Review measured and Li-informed organic curves, after 26a | Curve plots and tail assumptions |
| 26b | Update organic carbon, after 26a and 03/20 | Separate depth and curve effects |
| 26c | Assemble final organic + mineral total, after 26b and 22 | Component closure and mapping status |
| 27a | Support/extrapolation diagnostics, after 26c | Source and support maps; not confidence |
| 27b | Seven-scenario sensitivity rerun, after baseline | Sensitivity summaries; expensive |
| 28 | Individual-core holdouts, after baseline | Every held-out error and paired summaries |
| 29 | Optional conditional organic-carbon SD proxy, after 20/26a/26b | Proxy map and assumptions; not full posterior |
| 40, 41 | Community/technical report assets, after analysis and sensitivity | Asset manifests and narrative inputs |
| audit | Existing saved-output scientific checks, after complete workflow | Closure and source audits |
| record | Save hashes, configuration/code provenance and environment | Timestamped reproducibility folder |
| verify | Compare current files with latest record | Changed/missing/new files; does not rerun analysis |

Small tests are quick; raster mapping and auditing can take minutes or longer depending on the machine and disk. Sensitivity repeats large calculations. Check space on both output and temporary drives before 27b. Run `record` only after checking results; it records current files, it does not certify that all outputs are current. `verify` is an explicit provenance check, not a statistical validation test.

## Authoritative products

Under `outputs/route_l/`:

- `organic_horizon_thickness_blend.tif`: original and blended organic depth, source/influence diagnostics.
- `organic_carbon_li_community_update.tif`, band `updated_organic_carbon_kg_m2`: updated organic-only carbon; includes Li comparison and change decomposition bands.
- `full_variable_depth_carbon_components.tif`: final organic + provisional mineral total. Step 25's partial accounted output is not this total.
- `leave_one_core_out/`: individual-core holdouts. Different IDs do not establish independent sampling sites; no defensible grouped sampling units have been supplied.
- `organic_carbon_bayesian_uncertainty/`: **conditional SD proxy**, transferring the depth blend's fractional variance contraction to Li carbon spread. It omits carbon-curve uncertainty and is not a calibrated carbon posterior. See `docs/uncertainty_companion.md`.
- `simple_spatial_depth/`: earlier exploratory depth-model exports; not the production depth or carbon model, and not a prerequisite for Step 29.

Finite map coverage is not evidence of local accuracy; NoData is not confirmed absence of carbon, peat or water. The current geographic weight is exponential decay, `exp(-distance/range)`, rather than classical inverse-distance weighting.

## Documentation and checks

- `docs/run_in_RStudio.md`: practical running/review instructions.
- `docs/route_l_revision_review.md`: current baseline scientific rationale.
- `docs/uncertainty_companion.md`: optional proxy method and limits.
- `docs/candidate_covariates.md`: exploratory variables by carbon pool and landscape setting.
- `docs/alphaearth_local_pilot.md`: Colab/Python GEE covariate pilot using the uploaded AOI asset; no candidate fitted yet.
- `docs/covariate_analysis_locked_plan.md`: locked sequence and decision rules for the geographic, environmental, AlphaEarth, hydrology and external-model experiments.
- `config/covariate_analysis.R`: machine-readable covariate protocol; it does not alter the production Route L configuration.
- `tests/route_l_workflow_review.R`: portable runner and provenance checks.
- `tests/route_l_regression.R`, `tests/route_l_loo_review.R`, `tests/route_l_output_audit.R`: existing scientific checks.

Older planning/scaffold documents remain for history; they do not override current code, configuration or the documents above. S30/S100 are undeveloped scaffolds. Steps 40/41 regenerate figures, tables and Markdown, not the separately edited Word reports. Refresh those reports separately after analytical changes.

Before a covariate experiment, verify that the authoritative comparator is
unchanged with `Rscript scripts/50_freeze_covariate_baseline.R --verify`.
Validate the cloud raster exports and downloaded point tables with
`python3 scripts/gee/validate_covariate_exports.py`; results are written under
`outputs/covariate_analysis/part2_qa/`.

After downloading both the Riley/CanPEAT ecosystem batches and the community
core export, build reconciled analysis tables with
`Rscript scripts/51_prepare_covariate_point_tables.R`. Raw GEE CSVs remain
unchanged; prepared tables are written under `data/derived/covariate_analysis/`.
Run the outcome-blind full-vector distance and PCA diagnostic with
`Rscript scripts/52_diagnose_alphaearth_similarity.R`; results are written to
`outputs/covariate_analysis/part3_similarity_diagnostics/`. This step does not
use carbon responses or fit a carbon model.
Compare the locked geographic transfer with interpretable and AlphaEarth
similarity using whole-core holdouts with
`Rscript scripts/53_compare_local_similarity_holdouts.R`; results are written
to `outputs/covariate_analysis/part3_local_holdouts/`. This step does not write
a candidate raster or alter Route L products.
The local/cloud execution boundary and the conditions for using Vertex AI are
recorded in `docs/covariate_cloud_execution_plan.md`.
