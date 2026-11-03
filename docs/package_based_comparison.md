# Package-based carbon updates

Run `source("script_5_packages.R")` from the workshop project folder. The script reads the original inputs and existing handwritten predictions. It does not source, edit, or overwrite either handwritten script or its outputs. Results are written to `outputs/Testing_Bayes_Packages/stock_0_100/`.

## Package choices

| Analysis | Numerical engine | What it replaces |
| --- | --- | --- |
| Raw cores: gstat | `gstat::krige0()` | Handwritten Gaussian spatial conditioning and its variance calculation |
| Raw cores: geoR | `geoR::krige.bayes()` | A second package implementation of the same fixed-parameter spatial update |
| Area-based map | `mvtnorm::cond_mvnorm()` | Handwritten normal-normal precision weighting |
| IDW map | `gstat::idw()` and `mvtnorm::cond_mvnorm()` | Held-out IDW prediction and handwritten normal-normal precision weighting |

`geoR` 1.9-6 was installed in the R 4.6 library for this comparison; the other packages were already available. The run records package versions and full session information. `gstat` has peer-reviewed package papers; `mvtnorm` has published package literature and a Springer reference. The 2009 mvtnorm article discusses probability algorithms, not a review of this new workflow or specifically every current conditional-distribution function. The current conditioning interface is documented separately below.

`spBayes`, `brms`, `GauPro`, `GPfit` and `GPBoost` are useful alternatives for other spatial models. We did not select them as default fitted models here: changing the mean, covariance, observation error or parameter treatment would prevent a direct implementation comparison. `CARBayes` primarily addresses areal-neighbourhood models, which differs from the current point-location process. `fields::rdist()` can replace distance calculations, but is not by itself a replacement for the conditioning engine.

## Preserve the statistical model

The raw-core settings remain `range_m = 2000` and `observation_log_sd = 0.00`. No maximum-likelihood estimation, range fitting or new noise assumption was introduced.

The lognormal moment conversion remains:

```r
log_variance <- log1p((stock_sd / stock_mean)^2)
log_mean <- log(stock_mean) - log_variance / 2
```

The uncertainty raster is a Q95-Q05 prediction-interval width. Its division by `2 * qnorm(0.95)` is a working normal-interval approximation. Treating the resulting spatially varying variance as latent process uncertainty is also a model assumption. Source predictive uncertainty does not supply a separate decomposition into latent variance and measurement error.

Raw cores use the prior log-mean and log-variance at each location. Their residuals are standardised to a latent process with unit variance before calling `gstat`. Observation variances are divided by the corresponding prior variance. Predictions are transformed back using each target location's own prior mean and variance. This is algebraically equivalent to the handwritten heteroscedastic prior model.

This standardisation is important: the installed `gstat::krige0()` custom-covariance route uses the first target's variance when computing marginal variances. The standardised target process has variance 1 everywhere, so that assumption holds. Calling it directly with varying target variances would not preserve our model. We request marginal variances with `fullCovariance = FALSE`; no full posterior grid covariance or joint realisations are claimed.

### geoR comparison

`geoR::krige.bayes()` uses the same standardised raw-core residuals and prediction locations as `gstat`. Its settings are explicit:

```r
model = geoR::model.control(cov.model = "exponential", lambda = 1)
prior = geoR::prior.control(
  beta.prior = "fixed", beta = 0,
  sigmasq.prior = "fixed", sigmasq = 1,
  phi.prior = "fixed", phi = 2000,
  tausq.rel.prior = "fixed", tausq.rel = 0
)
```

The actual script takes `phi` and the relative nugget from the model settings rather than repeating constants. `lambda = 1` avoids transforming already logged and standardised data again. Output requests analytical predictive moments of the latent signal (`signal = TRUE`) and disables posterior and predictive simulations. Each pixel's prior log-mean and log-variance are then restored before lognormal back-transformation.

This is fixed-parameter Bayesian conditioning. It does not integrate uncertainty about the range, sill or nugget. The geoR documentation supports this configuration, while recommending `krige.conv()` for ordinary fixed-parameter prediction; we use `krige.bayes()` here to test the requested interface against the existing model.

geoR uses a constant relative nugget. Our common log-scale observation variance becomes `observation_log_sd^2 / core_log_variance` after standardisation. At the current zero observation error, these values are all zero and the two package models are equivalent. The script stops with an explanation if a future change produces heterogeneous standardised noise, rather than silently averaging it. The gstat implementation supports heterogeneous observation noise.

The geoR run is another implementation of the raw-core analysis, not an additional independent source of data. Area-based and IDW updates remain separate alternatives.

For each area-based or IDW pixel, `mvtnorm` conditions a bivariate Gaussian with standardised latent stock Z and noisy observation U:

```text
Var(U) = 1 + q, Var(Z) = 1, Cov(U, Z) = 1
q = working log-error variance / prior log-variance
```

The package computes the conditional mean and variance. The script performs input preparation, model specification, lognormal moment conversion, back-transformation, mapping and output comparison. Using packages reduces custom numerical code; those remaining modelling choices still require scientific assessment.

## Error assumptions

- Raw cores: fixed log-scale observation SD, currently zero for the synthetic test.
- Area-based map: between-core SD, about 19.05254 kg C/m², used as a working local-error scale. It is not standard error of the area mean or an independently calibrated map error.
- IDW: held-out-plot RMSE, about 9.305307 kg C/m², with Script 2's selected IDW settings fixed. This is working error calibration, not independent validation of the complete posterior update.
- Map-based likelihoods assume log-unbiased errors. The conversion `log1p((error_sd / mapped_estimate)^2)` is a plug-in approximation using the mapped value as reference.
- The three methods are alternatives derived from the same cores, not independent evidence sources to apply sequentially.
- Intervals are pointwise latent-stock intervals conditional on fixed assumptions. Area/IDW error covariance across pixels is not modelled, and these intervals do not provide regional-total uncertainty.
- Grid outputs are predictions at 1 km grid centres, not calculated 1 km block averages.

## Numerical comparison

The comparison uses all saved cells, requires identical cell sets and coordinates, checks grid geometry, and checks the raw-core baseline settings. It compares prior means/SDs, posterior means/SDs, both interval endpoints, changes from prior, and log-variance reduction. The tolerance is `1e-8 + 1e-8 * abs(handwritten_value)`.

| Method | Matched cells | Mean updated stock (kg C/m²) | Largest difference in posterior mean (kg C/m²) |
| --- | ---: | ---: | ---: |
| Raw cores: gstat | 12,578 | 18.18803 | 2.13e-14 |
| Raw cores: geoR | 12,578 | 18.18803 | 2.13e-14 |
| Area-based | 12,578 | 20.40827 | 2.49e-14 |
| IDW | 12,578 | 28.92105 | 4.97e-14 |

All 32 implementation/quantity comparisons passed. The direct geoR-versus-gstat check also passed for all six quantities; their largest posterior mean difference was 2.84e-14 kg C/m². The new geoR raster and CSV exports passed geometry, valid-cell count, finite-value, variance and interval-order checks. Differences are at floating-point rounding scale. This supports numerical equivalence of the implementations under the same assumptions; it does not demonstrate predictive accuracy or calibrated uncertainty.

Additional isolated checks confirmed raster export/reopening, interval ordering, the raw-core standardisation with heterogeneous nonzero observation errors, and mvtnorm conditioning across several signal/noise ratios. These checks did not change the published model settings.

## Outputs

Each of `raw_cores/`, `raw_cores_geoR/`, `area_based/` and `idw/` contains `predictions.csv`, `posterior_maps.tif`, `posterior_mean.png` and `differences_from_handwritten.csv`. The IDW folder also contains `error_checks.csv`.

The parent output folder contains:

- `comparison_summary.csv`: numerical comparison by implementation and quantity against the corresponding handwritten analysis. Both raw-core packages use the same handwritten raw-core reference.
- `geoR_vs_gstat.csv`: direct comparison of both raw-core package engines, including log-scale means and variances.
- `comparison_maps.png`: handwritten and package posterior means side by side.
- `difference_maps.png`: posterior mean differences, with values below 1e-8 displayed as zero; exact differences remain in the CSV files.
- `model_summary.csv`, `model_settings.csv`, `prepared_prior.tif`.
- `input_manifest.csv`, `package_versions.csv`, `package_references.txt`, `session_info.txt`.

For a new project without handwritten reference outputs, set `compare_handwritten <- FALSE`. This produces the package analyses without claiming a comparison.

## Notes on the additional background outline

The supplied four-phase outline is a useful teaching structure. However, Bayesian modelling does not generally require normal observations, and Gaussian-process models can use non-Gaussian observation likelihoods. Positivity and skewness motivate considering a log-Gaussian model, but logging data does not establish joint Gaussianity. The equations should preserve `SD^2 / mean^2`, `log_mean = log(mean) - log_variance / 2`, and the selected exponential kernel `exp(-distance / range)`; some pasted notation has lost fraction formatting.

Parameter estimation is an optional model extension, not a requirement for defining a conditional update. It remains deferred. Six widely spaced synthetic cores provide limited evidence for separating range and observation error.

## References and documentation

1. Pebesma, E.J. (2004). Multivariable geostatistics in S: the gstat package. *Computers & Geosciences*, 30, 683-691. https://doi.org/10.1016/j.cageo.2004.03.012
2. Gräler, B., Pebesma, E., and Heuvelink, G. (2016). Spatio-Temporal Interpolation using gstat. *The R Journal*, 8(1), 204-218. https://journal.r-project.org/articles/RJ-2016-014/
3. Mi, X., Miwa, T., and Hothorn, T. (2009). New Numerical Algorithm for Multivariate Normal Probabilities in Package mvtnorm. *The R Journal*, 1(1), 37-39. https://journal.r-project.org/articles/RJ-2009-001/
4. Genz, A., and Bretz, F. (2009). *Computation of Multivariate Normal and t Probabilities*. Springer. Use `citation("mvtnorm")` for the installed package's recommended citation.
5. gstat kriging and custom covariance API: https://r-spatial.github.io/gstat/reference/krige.html
6. mvtnorm conditional-normal API: https://search.r-project.org/CRAN/refmans/mvtnorm/html/margcond.html
7. Sothe et al. (2022), uncertainty definition, section 2.6: https://doi.org/10.1029/2021GB007213
8. Gramacy's recommended book citation is 2020: https://bookdown.org/rbg/surrogates/
9. Plant's supplied *Gaussian Processes in Spatial Data Analysis* supplement is dated 31 May 2025 and is stored in `docs/bayesian_spatial_analysis_resources/AdditionalTopics_GaussianProcesses.pdf`.

10. geoR Bayesian prediction and fixed-parameter controls: https://search.r-project.org/CRAN/refmans/geoR/html/krige.bayes.html
11. geoR output controls, analytical moments and latent-signal prediction: https://search.r-project.org/CRAN/refmans/geoR/html/output.control.html
12. Ribeiro Jr., P.J., and Diggle, P. (2025). *geoR: Analysis of Geostatistical Data*. Installed package version 1.9-6; recommended citation from `citation("geoR")`. https://doi.org/10.32614/CRAN.package.geoR
