# Workshop methods and fidelity audit

Audit of branch `workshop-mvp`, starting at `7b0d5a6`, against the original Word overview, earlier workshop materials, current Google explainer, maintained R code and bundled Carbon SOD summaries. Completed during the 26 October 2026 workshop revision.

## Original sequence and corrections

| Original component | Finding and correction |
| --- | --- |
| Compile field and laboratory measurements | Retained workbook and eight-column local table; entered OC remains authoritative. J7 = 0.09 g/cm³ is an import-only, identity-checked example correction. Outdated LOI, water and C:N chemistry remains flagged. No workbook bytes changed. |
| Stocks, depth profiles, SD and CV | Retained full slice/profile outputs, cumulative depth curves and complete-depth comparisons. Arithmetic and kg/m²-to-t/ha conversion checked independently. |
| Other ground datasets | Retained all five reference dimensions and seven useful figures. Summaries are aggregate records, not raw profiles; no hypothesis tests or causal claims are reconstructed. |
| Ecosystem and whole-area totals | Retained hypothetical nearest-core zones as a teaching example; added a real-data interface requiring explicit mapped ecosystem labels. Clipped polygon areas supply exact AOI area weights. Sampling representation remains an assumption. |
| Ordinary kriging | Missing in the previous implementation. Added an explicit, separate 240-point stationary simulation, training-only variogram fitting and buffered spatial assessment with `gstat`. Retained IDW as a simpler distinct method and replaced its production base-R interpolation with `gstat::idw`. |
| Raster/predictor analysis | Retained `stats::lm` and spatially held-out simulated observations; restored the original diagnostic orientation (predicted horizontal, observed vertical). Added local predictor/observation interfaces with explicit split and separation checks. |
| Compare methods and priors | Retained method totals, added a pre-update prior comparison on the same grid, and retained posterior comparisons. Separate simulations cannot rank the workbook maps’ real accuracy. |
| Bayesian updating | Kept the correct closed-form normal–normal calculation. The original draft equation adding means to variances is dimensionally invalid; precision weighting is required. SDs of 5 and 3 kg C/m² remain explicit teaching assumptions. |
| Case studies | Preserved author placeholders and restored the original planned prior/data-product list. Real Fort Severn analysis, original Li map specification, chronology and accumulation rates remain unresolved/deferred. |

## Scientific contracts

`stock = 0.1 × dry bulk density (g/cm³) × organic C (%) × thickness (cm)` gives kg C/m². Multiply by 10 for t C/ha. Integrate only measured slices; clip slices at reporting boundaries under a uniform-within-slice assumption. Complete depth coverage is required for 0–30 and 0–100 cm comparisons. Sampled totals over different depths are retained but not treated as comparable whole profiles. Local depth is measured from the ground surface including organic layers. Coarse-fragment and material-basis corrections must be resolved before import.

Reference concentration and density use thickness weights. Stock uses layer mass per area; these are different quantities. Reference prepared stocks may use different coarse-fragment conventions, estimated properties and a 0.5 cm coverage tolerance. Counts are profile records and may include repeated sites. Keep the limitations in `SOURCES.md`; these comparisons do not turn aggregates into raw observations or establish independent validation.

Area sums use exact clipped cell area. IDW and kriging values are point predictions represented on cells; they are not block means. Full GeoTIFF edge pixels extend outside the polygon. Boundary cell centres can lie outside the AOI; inventories approximate the continuous surface using those point values and exact clipped areas. Refining the grid tests numerical sensitivity but cannot create more observations. A community AOI may include water: use an appropriate land/carbon-pool mask for a real inventory.

AOI geometry is unioned and measured in an AOI-centred Lambert azimuthal equal-area projection. Distances use projected metres. This is useful for regional domains, but distances distort on continental extents; the code warns beyond a 1,000 km span. The national-extent test establishes software behaviour, not interpolation suitability. The 10,000-cell and two-million-distance-pair limits fail before allocating excessive matrices. Antimeridian-spanning inputs need explicit preprocessing.

The kriging simulation uses a known stationary isotropic exponential covariance, a constant mean of 18, nugget 1 and partial sill 9 (variance units). A separate training fit selects among five range candidates by variogram weighted residual error. `fit.variogram` fits nonnegative sills at fixed candidate ranges; singular fits are rejected. The holdout is the eastern quarter with a 5%-span buffer. The fixed run has 159 training, 60 holdout and 21 buffer observations. Kriging RMSE is 3.63 versus IDW 2.87 kg C/m²; nominal 95% intervals cover 80%. This is deliberately reported without seed selection to manufacture a favourable result. Conditional variance omits fitted-parameter uncertainty and does not validate the four-core map. Trend, anisotropy, alternative covariance forms, repeated spatial splits and adequate real sampling require scientific review.

For independent Gaussian errors, `Vpost = 1/(1/Vprior + 1/Vnew)` and `Meanpost = Vpost × (Meanprior/Vprior + Meannew/Vnew)`. Correlated errors, shared training data, bias and mismatched point/cell support invalidate a naive precision sum. Core-to-core variation is not calibrated map error. Cell SDs cannot be summed as if cell errors were independent to get total uncertainty. Normal support permits negative bounds; no silent truncation is applied.

## Sources consulted

- [gstat ordinary/point/block kriging and IDW reference](https://r-spatial.github.io/gstat/reference/krige.html): formula `z ~ 1` specifies ordinary kriging when the mean is not supplied; default prediction support follows the observations. The workshop deliberately uses point support.
- [gstat variogram-fitting reference](https://r-spatial.github.io/gstat/reference/fit.variogram.html): fitting can be singular and depend on starting values. The demonstration constrains ranges and checks fit status.
- [Pebesma 2004, Multivariable geostatistics in S](https://doi.org/10.1016/j.cageo.2004.03.012): methodological reference for the package.
- [Roberts et al. 2017, cross-validation for structured data](https://florianhartig.github.io/publications/roberts2017cross/): spatial assessment must match dependence and the prediction task. A single holdout here is a lesson, not comprehensive validation.
- [Natural Earth terms](https://www.naturalearthdata.com/about/terms-of-use/): public-domain cached context. The cache is independent of analytical land cover; source filenames and checksums are retained.

## Verification

The complete example was regenerated with R and gstat 2.1.6, sf 1.1.1, terra 1.9.27, readxl 1.5.0 and ggplot2 4.0.3. `outputs/R_session.txt` is generated on each complete run and is not versioned. Run `source("tests/run_checks.R")` after `source("run_all.R")` for reproducible checks.

Checks cover independently known stocks, unit conversion, zero carbon, invalid and missing measurements, overlaps, incomplete depth coverage, protected J7 correction, arbitrary ecosystem labels, gstat/base-arithmetic IDW agreement, duplicate locations, excluded/outside-AOI samples, western/eastern/northern Canadian AOIs, tiny and thin polygons, multipart and continental extents, irregular polygons with holes, projected boundaries, wrong longitude signs, grid/matrix limits, training/holdout separation, Bayesian limiting cases and failed reruns preserving learner notes. Local-interface tests exercise explicit missing-stage reporting, supplied ecosystem/predictor/prior data, old-grid rejection and failed-rerun cleanup. Synthetic test fixtures are removed after testing and never presented as observations.

The complete workflow now invalidates dependent saved outputs when an upstream script reruns. A failed script clears its partial products. Functions retain understandable arithmetic for stocks, summary statistics, linear regression and conjugate updating; gstat is the only added scientific dependency. Offline context adds no package dependency.

A fresh public branch ZIP downloaded with `usethis::use_course()` was unpacked under a temporary directory and ran `run_all.R` plus the entire check suite without the parent research project. The complete run retained its success marker. All 26 displayed Markdown code blocks parse, and each of the 20 figures matches its immediately following `view_figure()` path. The Google explainer has the same heading sequence and 20 inline figures; its exported PDF was visually checked, including keeping new kriging figures with captions and runnable code. The shorter own-data route also invalidates optional local results when the routes share an output directory.
