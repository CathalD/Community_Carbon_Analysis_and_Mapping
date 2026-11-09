# Route L covariate analysis: locked plan

**Protocol:** `route_l_covariate_experiment_v1`  
**Locked:** 21 September 2026  
**Clarified:** 23 September 2026 — external models are local-AOI prior/reference
learners only; DNN/stack candidates do not change the primary model order.  
**Baseline status:** authoritative and unchanged  
**Machine-readable decisions:** `config/covariate_analysis.R`

## Purpose

Test whether environmental information improves the transfer of community core
measurements across the AOI without silently replacing the existing Route L
maps. The existing geographic Route L outputs remain the comparator until a
candidate passes the stated adoption checks.

## Locked decisions

1. Use the complete 64-dimensional AlphaEarth vector for the first similarity
   experiment. L2-normalize nonzero vectors and use cosine/chord distance. PCA
   is a fallback only if distance-concentration or computing diagnostics show a
   real problem.
2. Retain the current conditional precision proxy. Do not add a naive effective
   sample-size or Voronoi adjustment. Continue to label the uncertainty output
   as a conditional SD proxy, not a calibrated posterior.
3. Compare geographic bandwidths of 500, 1,000, 2,000 and 4,000 m. The current
   2,000 m result remains the baseline. These are prespecified sensitivity
   cases, not a model-selection contest based on eight cores.
4. Leave the provisional mineral-carbon method unchanged and outside the first
   covariate experiment.
5. Build hydrology as a separate, auditable module. Begin with CDEM and
   Copernicus GLO-30 because currently checked public 1–2 m HRDEM/ArcticDEM
   products do not cover this AOI. Recheck Ontario and partner holdings before
   ruling out a local high-resolution DTM.
6. Do not impose a physical carbon-density clamp on accepted community data.
   Implied density may be reported as a diagnostic without changing outputs.
7. Train the first supervised external model on comparable Riley/CanPEAT
   complete 0–30 cm stocks. Keep the six comparable community profiles
   untouched for the final local transfer test. Do not call this a community-
   trained model. The external data inform a prior/reference relationship for
   the community AOI; they do not authorize a Canada-wide prediction map.
8. Defer anisotropic weighting until hydrology is defensible and residuals show
   directional structure. If tested, prefer terrain/flow-defined distance over
   fitting free anisotropy parameters to eight community sites.

## Work sequence

### Part 1 — Freeze the Route L baseline

- Record checksums for authoritative maps, QA reports, sensitivity and holdout
  summaries, raw community inputs and the current Route L configuration.
- Reference the earlier full input/code/output preservation record rather than
  repeatedly hashing the four multi-gigabyte Li source rasters.
- Verify the freeze before and after every candidate experiment.
- Never write candidate products over files in `outputs/route_l/`.

### Part 2 — Validate covariates

- Inventory downloaded raster groups and point tables.
- Check band names, units, native support, geometry, missingness and ranges.
- Confirm exact alignment to the 30 m Route L comparison grid.
- Compare raster values with exported community and external-core point values.

### Part 3 — Local similarity pilot

- Compare geographic-only, geographic × interpretable environment, and
  geographic × full-vector AlphaEarth weights.
- Run the four geographic bandwidth sensitivity cases.
- Diagnose AlphaEarth distance contrast before considering PCA.
- Report every held-out core alongside MAE, RMSE and bias.

### Part 4 — Supplemental hydrology and terrain

- Export candidate DEMs and quality masks over a hydrologically adequate buffer;
  do not rerun the full imagery stack.
- Use WhiteboxTools in small tutorial-style steps: inspect the DEM, condition
  depressions/flats, calculate a named flow direction and accumulation method,
  derive slope, specific catchment area and TWI, and inspect every intermediate.
- Compare CDEM and Copernicus results. Admit TWI only when flow paths and edge
  validity are credible.

### Part 5 — External 0–30 cm carbon model

- Train elastic net and random forest candidates on Riley/CanPEAT.
- Treat a small regularized DNN and a spatially valid stack as optional
  challengers, not required upgrades. The stack must use out-of-fold predictions
  from identical spatial/source-aware resamples.
- Use spatial blocks and study grouping to reduce leakage.
- Reserve all comparable community cores for the final transfer assessment.
- If transfer succeeds, test the regional prediction as a prior with a local
  community residual correction and predict only within the community AOI.

### Part 6 — Optional flow-aware transfer

- Consider hydrologic/cost distance only after the terrain module passes QA and
  residuals support an along-flow/across-flow distinction.
- Prefer hydrological covariates alone if they capture the useful structure;
  explicit anisotropy is not required.

### Part 7 — Decision and reporting

- Adopt a candidate only if improvements are meaningful, per-core failures are
  visible, results are reasonably stable, coverage/support is acceptable and
  the map has a defensible landscape interpretation.
- Otherwise retain Route L and use the results to identify environmental and
  geographic sampling gaps.
- Revisit uncertainty and reports only after a candidate passes this gate.

## WhiteboxTools learning approach

The hydrology work will be taught as a reproducible sequence rather than one
wrapper function. Each operation will document its input, output, parameters,
physical meaning, common failure modes and a visual or numerical check. The
workflow will retain intermediate rasters so that depression treatment, flat
resolution and flow routing can be understood and challenged independently.

## Current checkpoint

**Part 1 completed on 21 September 2026.** The targeted freeze contains 26
authoritative files and references the earlier full SHA-256 preservation
record. The independent verification command passed. Headline baseline values
and limitations are recorded beside the manifest in `SUMMARY.md`.

**Part 2 completed on 23 September 2026.** All 56 automated checks passed.
The audit covers 12 cloud COGs (21.58 GiB), 119 predictor bands, five QA bands,
the eight community rows and 1,564 Riley/CanPEAT rows. All 992 community
point/band values matched the native cloud pixels, including paired missing
values. The largest absolute finite difference was 2.47e-7.

**Part 3 completed on 24 September 2026.** The outcome-blind 64-D AlphaEarth diagnostic found
substantial distance contrast: among 1,009 unique valid external locations,
the median chord distance was 0.906 and the 5th–95th percentile range was
0.436–1.304. PCs 1–5 explained only 71.9% of variance, so the evidence does
not justify replacing the full vector with 3–5 PCs. The primary environmental
bandwidth is frozen outcome-blind at the unique-location median nearest-
neighbour distance (0.2325). All 12 geographic-only 2,000 m predictions
reproduced the existing Route L holdouts within 0.0001.

Neither environmental-similarity candidate improved the locked 2,000 m
comparator. Thickness MAE was 15.48 cm for geographic-only, 28.57 cm with
AlphaEarth and 19.42 cm with the interpretable metric. Organic-carbon MAE was
1.237, 1.310 and 1.388 kg C/m2, respectively. AlphaEarth also fell below the
0.1 equivalent-core support threshold for FS-01 and FS-02 even at 4,000 m.
Geographic-only 4,000 m reduced thickness MAE to 9.27 cm while organic MAE was
essentially unchanged at 1.214 kg C/m2, but this is retained as a sensitivity
finding rather than selected from eight/four holdouts. No environmental
candidate advances to full-AOI mapping; Route L remains authoritative.

**Part 4 setup completed on 24 September 2026.** The standalone GEE exporter
`scripts/gee/export_hydrology_dems.py` prepares buffered CDEM, Copernicus
GLO-30 and Copernicus QA inputs without rerunning the imagery stack. The local
runner `scripts/54_build_whitebox_hydrology.R` has explicit, retained stages
for least-cost breaching, flat/depression resolution, slope, D-infinity flow
direction, unlogged specific catchment area and strict TWI. Processing is
waiting only for the two DEM files to be exported from Cloud Storage and placed
under `data/raw/gee_imports/hydrology/`; no hydrology result has yet passed QA.

**Part 4 processing and audit completed on 25 September 2026.** Both DEMs were
processed successfully after conversion from Earth Engine COGs to value-checked
Float32 GeoTIFFs. Coverage is effectively complete and the Copernicus quality
masks are internally consistent. Conditioned elevation agrees strongly between
DEMs (r=0.998), but specific catchment area (r=0.157) and strict TWI (r=0.679)
are much less stable across the AOI. At the eight community cores, TWI is
negatively correlated between candidates (r=-0.362), with a median absolute
difference of 1.80 and maximum of 6.96. Raw TWI therefore does not pass the
admission gate and will not enter the carbon model at this checkpoint. Detailed
results are recorded in `docs/hydrology_part4_findings.md`.

**Part 5 completed on 25 September 2026.** Elastic net and random forest models
were trained on 1,512 complete Riley/CanPEAT 0–30 cm profiles using 105
community-complete covariates. Five-fold nested validation kept shared 100 km
blocks, source studies and coordinates together, and coordinate weights
prevented repeated profiles from dominating. At 998 held-out unique locations,
the random forest improved MAE by 5.7% over the external-mean comparator (3.598
versus 3.816 kg C/m2), but explained only 10.4% of held-out variance. Elastic
net did not improve the comparator.

The untouched six-profile community transfer test failed. Random-forest MAE
was 5.779 kg C/m2, versus 4.530 for the already biased external mean, and its
bias was +5.276 kg C/m2. Therefore no external model advances to local-prior
blending, cloud deployment or AOI raster prediction. Part 6 remains deferred;
the workflow proceeds to the Part 7 decision and reporting checkpoint. Full
methods and per-core evidence are in `docs/external_model_part5_findings.md`.
