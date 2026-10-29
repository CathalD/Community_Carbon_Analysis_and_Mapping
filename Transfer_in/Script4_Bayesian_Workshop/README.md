# Part 4 — Updating a carbon map with core observations

This workshop continues Scripts 1 and 2. We take the carbon stocks already calculated for each core, read an existing carbon map and its uncertainty, and update the map using a Bayesian spatial model.

The script follows the same sequence throughout: read the data, inspect it, run the calculation, make maps, and save the results. The statistical calculation is kept in one small function so that it can be explained separately.

**Use [Script_4.R](Script_4.R) with the existing `Transfer_Oct` example.** It does not generate new cores, carbon values or prior maps. The supplied `AOI-SIM-*` records already carry a simulated teaching-data label in the source files; that provenance is preserved in the exported core tables. These are the existing workshop records, not newly generated example data.

## What we mean by a depth-specific stock

We use the existing Script 1 columns without recalculating them:

| Selected depth | Input column |
|---|---|
| 0–30 cm | `C Stock 0-30cm` |
| 0–50 cm | `C Stock 0-50cm` |
| 0–100 cm | `C Stock 0-100cm` |

Script 1 sums section stocks whose **bottom depth is at or before the selected maximum depth**. Sections deeper than that maximum are excluded. Under the workshop assumption that each core represents its full profile, a profile ending at 60 cm still contributes its summed stock as a 0–100 cm observation. We do not remove it for being shallower than 100 cm, extend it to 100 cm, or include sections deeper than 100 cm.

For example, the supplied `AOI-SIM-C01` has a reported maximum depth of 30 cm and a `C Stock 0-100cm` value of 12.33 kg C/m². Script 4 retains 12.33 when mapping 0–100 cm. This explicitly follows the workshop's full-profile assumption, even where older notes in the source table describe incomplete recovery.

A section crossing a cutoff is handled exactly as it was in Script 1: a section ending below the cutoff is excluded rather than apportioned. Any future change to that calculation belongs in Script 1, so that all mapping scripts use the same stocks.

## Before running

Place this folder inside `Transfer_Oct` and use `Transfer_Oct` as the R working directory:

```text
Transfer_Oct/
├── Script_1.R
├── Script_2.R
├── aoi.geojson
├── Draft1/data/mapping_data.csv
├── data/priors/
│   ├── carbon_0_30_mean.tif
│   └── carbon_0_30_uncertainty.tif
└── Script4_Bayesian_Workshop/
    ├── Script_4.R
    └── README.md
```

The two TIFF filenames are the expected locations for exports of **your existing GEE assets**. Update the paths if your exports have different names. Export the mean and uncertainty bands for the same depth, carbon pool, release, extent and pixel grid, one band per TIFF. Keep NoData masked. Use a mean stock map, not carbon concentration, peat thickness, or a median mislabeled as a mean.

The maps should represent the same surface-based depth window and carbon pool as the core stocks. A peat-only stock and a total soil stock that includes additional mineral-soil carbon are different quantities.

Install the four packages once:

```r
install.packages(c("sf", "terra", "ggplot2", "patchwork"))
```

Edit Step 1, then run the script a section at a time in RStudio, or run:

```r
source("Script4_Bayesian_Workshop/Script_4.R")
```

## Step 1 — Choose the depth, inputs and assumptions

Start with `target_depth_cm <- 30`. This selects the corresponding Script 1 column. Set the two prior paths and `prior_source`, including the publication, release, GEE asset identifier and uncertainty definition.

| Setting | Meaning |
|---|---|
| `uncertainty_type` | Whether the uncertainty raster holds SD, variance, or a full prediction-interval width |
| `mean_to_kg_m2` | Conversion from source mean units to kg C/m² |
| `uncertainty_to_kg_m2` | Linear unit conversion for SD or interval width; for variance, applied after taking its square root |
| `prior_sd_multiplier` | Multiplier for a sensitivity run; 1 leaves the converted SD unchanged |
| `range_m` | Distance at which spatial correlation falls to about 0.37 |
| `observation_log_sd` | Assumed observation and local representativeness SD on the natural-log scale |
| `cell_size_m` | Prediction-grid spacing; initially 1,000 m, matching the scale used in Script 2 |
| `output_folder` | Destination for this depth and set of assumptions |

The initial 5 km range and 0.15 observation SD are **workshop assumptions**, not fitted estimates or values established by the cited papers. A log SD of 0.15 represents a multiplicative one-SD factor of about 1.16. It allows for differences between a core and the mapped field as well as measurement variation; the stock values themselves remain unchanged.

For t C/ha, use a conversion of 0.1 to kg C/m². For kg C/m², use 1. Rerunning in the same output folder replaces the named outputs. Use separate folders when comparing assumptions.

## Step 2 — Load the packages

| Package | Role in this script | Documentation and citation |
|---|---|---|
| `sf` | Read the AOI, create core points, transform coordinates and check containment | [sf documentation and citation guidance](https://r-spatial.github.io/sf/) |
| `terra` | Read, crop, sample, project and write rasters | [CRAN package page](https://cran.r-project.org/package=terra), DOI [10.32614/CRAN.package.terra](https://doi.org/10.32614/CRAN.package.terra) |
| `ggplot2` | Draw the core map and raster figures | [Official documentation](https://ggplot2.tidyverse.org/) |
| `patchwork` | Combine four maps into one figure | [Official documentation](https://patchwork.data-imaginist.com/) |
| R `base` and `stats` | Distances, matrix calculations, Normal quantiles and summaries | [Cholesky calculation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/chol.html), [lognormal distribution](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Lognormal.html) |

The Bayesian calculation uses R's matrix operations. There is no separate fitting package to install for this conditional Gaussian model. The script saves `citation()` results for the installed packages and a `sessionInfo()` record. Coordinate transformations use installed PROJ resources without downloading additional grids during the workshop.

## Step 3 — Read the same core table

`read.csv(..., check.names = FALSE)` preserves the column names exported by Script 1. We select the requested stock column and give it the shorter name `stock`.

The table retains `Core ID`, `Plot ID`, longitude, latitude and source notes. `Plot ID` identifies observations that should be held out together during validation. The workbook is not reread and stock values are not recalculated.

**Check:** look at the printed core IDs and stocks. Changing the target to 100 should select the existing 0–100 cm values, including stocks from shallower full profiles.

## Step 4 — Check the locations

The AOI and core points are transformed to EPSG:3978, a Canadian projected coordinate system with metre units. This lets us calculate distances in metres.

The script checks IDs and removes records with invalid coordinates, nonpositive or missing stocks, or locations outside the AOI. A later check also removes cores without usable prior values. Reasons are written to `excluded_cores.csv`. The positive-stock condition is necessary because this model takes logarithms; genuine zero-carbon observations would require a different likelihood.

There is no minimum core-depth filter.

## Step 5 — Read the prior and interpret its uncertainty

The script checks that the mean and uncertainty rasters share a grid and CRS, crops them to the AOI bounding box, then converts both to the required units. At each core, it samples the native prior cell before creating the output grid.

The Bayesian model needs an SD. The accepted conversions are:

| Source uncertainty | `uncertainty_type` | Conversion to SD before any unit multiplier |
|---|---|---|
| Standard deviation | `"sd"` | Use directly |
| Variance | `"variance"` | Square root |
| Full 90% interval width, upper minus lower | `"pi90_width"` | Divide by `2 * qnorm(0.95)`, approximately 3.29 |
| Full 95% interval width | `"pi95_width"` | Divide by `2 * qnorm(0.975)`, approximately 3.92 |

The width conversions are **Normal-equivalent approximations**. They do not recover an arbitrary skewed distribution from two quantiles. Do not use them for an interval half-width, coefficient of variation, standard error, or relative uncertainty index without adapting the interpretation.

For the existing Sothe/McMaster/WWF Canada 0–1 m product, the paper describes uncertainty derived from Q95–Q05 ranges, including SOC and bulk-density uncertainty. It should not be treated directly as one SD. This workshop uses `pi90_width` as an explicit approximation for that product. See [Sothe et al. (2022), uncertainty methods](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021GB007213).

We treat the converted map uncertainty as uncertainty about the latent stock field. Published prediction uncertainty may also contain small-scale and observation variation, so this allocation is a modeling assumption. Do not count the same error again in `observation_log_sd`. Ensemble disagreement alone is not automatically a calibrated stock SD.

## Step 6 — Make the prediction grid

The AOI is covered with regularly spaced prediction locations. Nearest-cell resampling keeps each selected prior mean paired with its uncertainty. Cells without valid positive prior means and SDs remain NoData.

These are **point predictions at grid centres**, using the native raster values as local prior summaries. They are not averages over each 1 km square. Changing resolution does not create more field information. Explicitly reconciling core support with raster-area averages is an extension discussed below.

## Step 7 — Put the prior on the log scale

Carbon stocks are positive. We therefore describe log stock with a Gaussian distribution, giving a lognormal stock distribution.

If a prior cell has arithmetic mean `M` and SD `S`, we set:

```r
log_variance <- log1p((S / M)^2)
log_mean <- log(M) - log_variance / 2
```

These formulas preserve the supplied mean and SD under the lognormal assumption. They follow from the [lognormal moment formulas in R's documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Lognormal.html).

The observation model is `log(observed stock) = log(true stock) + Normal(0, observation_log_sd^2)`. Errors are centred on zero in log units. On the stock scale this makes observations median-unbiased, rather than arithmetic-mean-unbiased.

## Step 8 — Update the spatial distribution

The function `condition_spatial()` performs the Bayesian update. It uses the discrepancy between the observed core stock and its prior expectation to adjust surrounding locations. The assumed spatial correlation is:

```r
correlation <- exp(-distance / range_m)
```

Nearby observations share information. Far from the cores, the map tends back toward its prior distribution. The calculation accounts for correlation between cores, rather than treating nearby cores as independent information about every pixel.

<details>
<summary>The equations behind the function</summary>

Let `z` be observed log stocks, `m` their prior log means, and `v` their prior log variances. At core locations the spatial covariance is

$$K_{ij}=\sqrt{v_i v_j}\exp(-d_{ij}/L).$$

Write `R` for the diagonal observation-error covariance, with entries `observation_log_sd^2`, and `k*` for covariances between a prediction location and the cores. The conditional log-stock mean and variance are

$$m_*^{post}=m_*+k_*^T(K+R)^{-1}(z-m),$$

$$v_*^{post}=v_*-k_*^T(K+R)^{-1}k_*.$$

The code uses `chol()`, `forwardsolve()` and `backsolve()` to evaluate these equations without explicitly inverting a matrix. This is the Gaussian-process conditioning calculation described by [Rasmussen and Williams (2006), Chapter 2](https://gaussianprocess.org/gpml/chapters/RW2.pdf), with a spatially varying prior mean and variance supplied by the maps.

</details>

The map SDs provide marginal uncertainty, but do not identify spatial correlation. The exponential correlation and its range are additional assumptions. The update is exact **conditional on these chosen assumptions**; it does not estimate their posterior distributions or fit an AOI-wide bias term.

## Step 9 — Predict stock and uncertainty

Predictions are calculated in batches to avoid holding every grid-to-core distance in memory. Each conditional log distribution is converted back to stock units:

```r
posterior_mean <- exp(posterior_log_mean + posterior_log_var / 2)
posterior_median <- exp(posterior_log_mean)
```

The script also calculates posterior SD and the 2.5th and 97.5th percentiles. These are intervals for the latent stock at each prediction location, conditional on the model assumptions. Exponentiating the log mean alone gives the median, not the arithmetic mean.

Log variance cannot increase in this fixed Gaussian update. SD in kg C/m² can increase if the estimated stock rises sufficiently, so a higher posterior SD is not automatically an error.

## Step 10 — Leave a sampling group out

For each `Plot ID`, withhold its cores, update using the remaining plots, and predict the withheld stocks. The output compares prior and updated predictive means using MAE, RMSE and bias, and reports the fraction of observations inside the updated 95% predictive intervals.

Here the prediction is for a **new observed core**, so we add observation variance to the latent field variance. Both prior and updated means use that same observation model. They therefore differ slightly from the latent mean maps.

The prior must be independent of these validation observations. A Script 2 map fitted to the same cores cannot serve as an independent prior. Leaving cores out of this update does not remove them from an earlier map's training set. With six teaching records, the scores illustrate the procedure; they do not establish performance on new field data. Plot holdouts also do not establish accuracy far beyond sampled locations.

## Step 11 — Save the maps

`bayesian_carbon_maps.tif` contains these bands:

| Bands | Interpretation |
|---|---|
| `prior_mean`, `prior_sd` | Converted prior summaries sampled on the output grid |
| `posterior_mean`, `posterior_median`, `posterior_sd` | Updated latent stock summaries |
| `lower95`, `upper95`, `interval_width95` | Conditional 95% stock interval and its width |
| `change_from_prior` | Posterior mean minus prior mean |
| `log_variance_reduction_fraction` | Fraction of prior log variance removed by conditioning |
| `nearest_core_m` | Distance to the nearest included core, in metres |

All stock quantities use kg C/m². Variance reduction is dimensionless. NoData marks cells that were not predicted.

This script does not calculate an inventory total or total uncertainty. Marginal pixel SDs are insufficient for an AOI-total interval because predictions are spatially correlated; a later inventory calculation needs joint covariance or joint draws and appropriate area support.

## Step 12 — Compare four maps

The saved figure shows the prior mean, updated mean, change in mean and updated interval width. The two mean maps share a colour scale. White points show the included cores.

**Discuss:** where did the map change most? How does that relate to core disagreement, prior uncertainty and distance? A reduction in conditional variance is a model consequence, not proof that the chosen prior or observations are accurate.

## Step 13 — Save a record of the run

The output folder also contains `included_cores.csv`, `excluded_cores.csv`, `heldout_predictions.csv`, `validation_scores.csv`, an RDS results file, settings, input checksums, package citations and session information. These record the stock column, full-profile assumption and uncertainty choices alongside the maps.

## How this uses the textbook

The supplied [Chapter 9 of Moraga's *Spatial Statistics for Data Science*](https://www.paulamoraga.com/book-spatial/bayesian-spatial-models.html) introduces hierarchical spatial models for areal observations, including the BYM model. We carry forward its separation of observations, an underlying spatial process and prior assumptions.

For cores at coordinates, [Chapter 15, model-based geostatistics](https://www.paulamoraga.com/book-spatial/sec-geostatisticaldataSPDE.html), is the closer match. Our model is `log(stock) = prior log mean + spatial deviation + observation error`. The map supplies the prior mean and marginal uncertainty; a Gaussian spatial field links locations.

This workshop computes the conditional Gaussian update directly. A fuller INLA–SPDE implementation would represent the field on a mesh and specify priors for spatial parameters. Point and area observations would need suitable projection operators. That would extend the model beyond its fixed range and point-support approximation; simply replacing a function call with `inla()` would not do so.

## Research connecting maps and ground observations

These studies motivate the approach; the script is not a reproduction of their complete models.

| Study | What was combined and how it relates |
|---|---|
| [Fuentes & Raftery (2005), *Model Evaluation and Spatial Interpolation by Bayesian Combination of Observations with Outputs from Numerical Models*](https://sites.stat.washington.edu/raftery/Research/PDF/fuentes2005.pdf). Biometrics 61, 36–45. DOI: 10.1111/j.0006-341X.2005.030821.x | Ground pollution observations and numerical-model output informed a latent spatial field, with model bias and differing spatial resolutions addressed. A foundation for Bayesian map–observation fusion. |
| [Verdin et al. (2015), *A Bayesian kriging approach for blending satellite and ground precipitation observations*](https://civil.colorado.edu/~verdina/Verdin-et-al-WRR-2015.pdf). Water Resources Research 51, 908–921. DOI: 10.1002/2014WR015963 | Gauges were related to satellite estimates and elevation through a Bayesian model, followed by spatial prediction. This supports map-informed prediction with ground-based corrections; it does not simply assimilate a published uncertainty raster as we do here. |
| [Moraga et al. (2017), *A geostatistical model for combined analysis of point-level and area-level data using INLA and SPDE*](https://researchportal.bath.ac.uk/en/publications/a-geostatistical-model-for-combined-analysis-of-point-level-and-a/). Spatial Statistics 21, 27–41. DOI: 10.1016/j.spasta.2017.04.006 | Point and areal data were linked to a shared field using support-aware projections. This is the main extension when raster cells must be modeled as area averages. |
| [Pallegedara Dewage, Minasny & Malone (2020), *Disaggregating a regional-extent digital soil map using Bayesian area-to-point regression kriging for farm-scale soil carbon assessment*](https://soil.copernicus.org/articles/6/359/2020/). SOIL 6, 359–369. DOI: 10.5194/soil-6-359-2020 | A soil-carbon example of Bayesian downscaling of an existing map. Its ground measurements were used for validation; it should not be described as the same ground-data updating model implemented here. |

## Repeat for another depth or assumption

For 0–100 cm, set `target_depth_cm <- 100`, choose the matching mean **and** uncertainty rasters, and change the output folder. All six supplied core records remain eligible under the full-profile assumption, subject to the spatial and positive-value checks.

Run each depth separately. Do not substitute 0–100 cm uncertainty for missing 0–30 cm uncertainty, combine several correlated prior maps as independent evidence, or add cumulative 0–30 and 0–100 cm stocks together. Separate fits do not enforce increasing stock with depth or estimate cross-depth covariance.

For a short workshop exercise, repeat the same depth with `range_m` values of 2,500, 5,000 and 10,000. Then try a larger `prior_sd_multiplier` or observation SD while keeping everything else fixed. Compare changes in the maps and held-out results. Selecting the best setting on these same holdouts would require separate evaluation before claiming an accuracy improvement.

## Verification with the supplied files

The script was checked with the existing six-core `mapping_data.csv`, the supplied AOI, and the locally available McMaster/WWF Canada version 3.0 **0–100 cm** mean and uncertainty maps, using the interval-width approximation described above. The default 0–30 cm run still requires the matching uncertainty export from your GEE assets; an uncertainty map was not invented for it.

The checks passed: all six supplied 0–100 cm stocks were retained, maps and held-out outputs were produced, interval bounds were ordered, conditional log variances were nonnegative and no greater than their priors, and the conditioning function agreed with an independent matrix calculation. The four-panel figure was also inspected.

For that run, prior and updated held-out RMSE were approximately 26.53 and 26.40 kg C/m²; only two of six observations fell inside the updated 95% predictive intervals. This is a working calculation, not a well-calibrated final mapping model. The supplied teaching observations disagree substantially with this national prior, and the assumed 5 km correlation provides little information between widely separated held-out cores. The test verifies the workflow, not ecological validity or improved field performance.
