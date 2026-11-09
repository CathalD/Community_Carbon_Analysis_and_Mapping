# Carbon Analysis and Mapping Workshop

Open `Carbon_Analysis_and_Mapping_Workshop.Rproj` in RStudio and run the numbered
scripts in order. Paths are relative to the project folder. Shared paths, CRS,
cell size and depth settings are in `config.R`; method settings remain near the
start of the relevant script.

## Inputs

- Field workbook: `data/field/field_cores.xlsm`.
- Study boundary: `data/boundary/aoi.geojson`.
- Reference observations and source notes: `data/reference/`.
- Script 3 covariates: `data/covariates/topography.tif` and `water.tif`.
- Prior stock and Q95-Q05 uncertainty width: `data/priors/`.

The supplied six-core workbook and AOI are simulated teaching data. Existing
results from before the cleanup remain on disk; rerun the scripts to generate
the products described below.

## Run the workshop

| Script | Purpose | Current outputs |
| --- | --- | --- |
| `script_1.R` | Field stocks, profile plots and reference context | `outputs/script1/`, including the mapping table |
| `script_2.R` | Shared grid, area-average log model, IDW and ordinary kriging | `data/grid/` and prediction, variance and log-moment rasters under `outputs/script2/` |
| `script_3.R` | Existing covariate analysis and model selection | Existing native-resolution products under `outputs/script3/<stock variable>/` |
| `script_4.R` | Prepare the prior once and compare matching-depth maps | Aligned prior moments, regression copy and comparison rasters under `outputs/script4/` |
| `script_5.R` | Separate raw-core, area-based and IDW Bayesian updates | Prediction and variance rasters under `outputs/script5/stock_0_<depth>/<method>/` |

Each script prints its results, parameter estimates and relevant warnings to the
console. Users can change the settings and rerun. No extra QA-report workflow is
required. Field inputs and generated outputs should be kept separate: update the
workbook, rerun Script 1, then rerun the affected later scripts.

## Depth convention

Script 1's stock calculations are unchanged. Every supplied core contributes to
every depth bin: integrate to the smaller of the target depth and recorded core
maximum. The field profile is assumed complete at its recorded maximum, as
specified for this workshop. No exclusion based on core maximum is applied in
the mapping or Bayesian scripts. The example therefore retains six cores in
0-30, 0-50, 0-100 and sampled-profile total groups.

The external reference collections retain their own documented depth eligibility
rules. Their source tables are not rewritten by this cleanup.

## Shared grid

Script 2 creates one template anchored to whole multiples of `cell_size_m` in
`analysis_crs`, with a saved boundary mask and stable cell IDs. All final spatial
interpolation and Bayesian products use this geometry and the same mask rule:
retain cells intersecting the boundary. `data/grid/grid_settings.csv` records the
boundary checksum and settings; `session_info.txt` records the software versions.
Changing the boundary, cell size or CRS requires rerunning the dependent scripts.

Script 3 deliberately keeps its original native-grid terrain calculations,
predictors, formulas, response scale, validation, model selection and mapping.
Script 4 transfers its selected map to the shared grid for comparison. External
prior stock and uncertainty are transferred together using nearest source values
at target cell centres. This is not uncertainty aggregation for block averages.

## Prediction scales and uncertainty

Script 2 fits its interpolation models to natural-log carbon. Power zero gives
the equal-weight mean of log C; positive powers give distance-weighted log C.
This changes the old arithmetic-stock baseline and raw-stock IDW analysis.
Stocks are strictly positive after exponentiation, but logging alone does not
establish Gaussianity or a well-supported variogram.

For area average and IDW, whole plots are withheld to estimate log errors. IDW
settings are selected within each outer training fold. The mean log error is
removed from the mapped log location, and the centered log-error variance is
used as a working likelihood variance. These estimated calibration parameters
are treated as fixed; the small synthetic example does not establish calibrated
uncertainty for a real landscape.

The exported prediction is the working lognormal mean `exp(mu + v/2)`; variance
is `expm1(v) * exp(2*mu + v)`. It is not a geometric mean or a standard error of
the landscape mean. Kriging uses its conditional log-scale prediction moments
with the same lognormal approximation. Fit warnings remain visible in the
console; a returned curve or raster is not proof that the fit converged.

For each method and depth, Script 2 exports:

- `<method>_<depth>cm.tif`: mean stock, kg C/m2.
- `<method>_variance_<depth>cm.tif`: pointwise variance, (kg C/m2)^2.
- `<method>_log_moments_<depth>cm.tif`: log location and log variance for Script 5.
- Area/IDW calibration and held-out tables: the inputs behind the error estimate.

The prior's Q95-Q05 interval width is approximated as stock SD by dividing by
`2*qnorm(0.95)`, then moment-matched to a lognormal distribution in Script 4.
The source predictive uncertainty is used as prior process uncertainty as a
working model assumption. Script 5 does not convert it or resample it again.

Script 5 preserves the raw-core exponential covariance model and its fixed range
and observation-error settings. Area and IDW updates read Script 2's calibrated
log moments directly. Each method updates the original prior separately; these
are alternatives based on the same field evidence. Pointwise variances do not
provide uncertainty for a study-area total without cross-cell covariance.

Script 3 is the explicit exception to the log-model change: its original analysis
and response scale are retained by request.

## Optional development material

`Testing_Bayesian_Models.R` forwards to the main Script 5, so there is only one
handwritten implementation to maintain. `script_5_packages.R` optionally compares
the current numerical update with gstat and mvtnorm and prints differences.
Original development scripts are retained under `archive/development/`.

Original transfer material remains under `archive/transfer_original/`. Source
terms, attribution and unresolved provenance are documented with the reference
data. No third-party data licence is changed by this cleanup.
