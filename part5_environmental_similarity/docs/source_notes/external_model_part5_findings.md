# Part 5 external 0–30 cm carbon model findings

**Run date:** 25 September 2026  
**Protocol:** `route_l_external_model_v1`  
**Decision:** do not use the fitted external model as a local prior or map it

## Question

Can a covariate model trained only on Riley and CanPEAT complete 0–30 cm carbon
stocks provide a useful regional prior for the six comparable community cores?
This was a transfer test, not a Canada-wide mapping exercise and not a model
trained on the community outcomes.

## Data and predictors

- 1,512 external profiles: 1,211 CanPEAT and 301 Riley.
- 998 unique coordinates and 887 source studies.
- Six community profiles with complete measured 0–30 cm stocks were withheld
  from all fitting and tuning.
- The primary feature set contains 105 covariates that are available at all six
  community profiles: 64 AlphaEarth dimensions, 35 surface-context bands and
  six terrain bands.
- Fourteen spring NDMI/TCW bands were excluded because one complete community
  profile lacks them. Raw Whitebox TWI was not admitted after the Part 4 DEM
  sensitivity audit.
- Missing external predictor values were filled using training-fold medians;
  centring, scaling and zero-variance removal were also learned separately
  inside each training fold.

## Leakage controls

Records sharing a 100 km projected block, source study or rounded coordinate
were joined into one connected validation group. The 180 resulting groups were
assigned to five balanced outer folds. Every outer training set was tuned again
with grouped inner folds. Repeated profiles at one coordinate received weights
that sum to one, so duplicate or layered records could not dominate training or
the primary coordinate-balanced evaluation.

All checks passed: a coordinate, study, spatial block and connected validation
group occurs in exactly one outer fold.

## Models

- Coordinate-weighted external training mean, as the null comparator.
- Elastic net (`glmnet`), tuning alpha and lambda.
- Random forest (`ranger`), tuning `mtry`, minimum node size and sample fraction.

The final external-only selections were alpha 0.25 and lambda 2.5929 for the
elastic net, and `mtry = 10`, minimum node size 15 and sample fraction 0.8 for
the random forest.

## External grouped holdout results

The most relevant summary averages repeated profiles at each held-out location.

| Model | MAE (kg C/m2) | RMSE (kg C/m2) | R-squared |
|---|---:|---:|---:|
| External training mean | 3.816 | 5.573 | -0.005 |
| Elastic net | 3.843 | 5.662 | -0.038 |
| Random forest | 3.598 | 5.262 | 0.104 |

The random forest therefore found a real but weak external relationship: its
MAE and RMSE were about 5.7% and 5.6% lower than the external-mean comparator.
The elastic net did not improve the null model. The random forest's low
R-squared shows that most held-out variation remained unexplained.

## Untouched community transfer test

| Model | MAE (kg C/m2) | RMSE (kg C/m2) | Bias (kg C/m2) | R-squared |
|---|---:|---:|---:|---:|
| External training mean | 4.530 | 5.072 | +3.975 | -1.592 |
| Elastic net | 4.612 | 5.111 | +4.057 | -1.633 |
| Random forest | 5.779 | 6.907 | +5.276 | -3.808 |

Observed community stocks range from 4.22 to 14.26 kg C/m2 and average 8.62.
The coordinate-weighted external mean is 12.60 kg C/m2. The elastic-net
predictions are almost constant (12.24–13.43), while the random forest predicts
12.17–17.09 and badly overpredicts the two lowest community stocks. Relative to
the already poor external-mean transfer, the random forest increases community
MAE by 27.6% and RMSE by 36.2%.

This is evidence of a target/population transfer mismatch. The covariates carry
some broad external information, but the relationship learned from Riley and
CanPEAT is not calibrated for the local 0–30 cm community stocks.

## Decision

Part 5 fails the local-transfer adoption gate. Neither fitted model will be
used as a regional prior, blended with Route L, sent to Vertex AI, or applied
to AOI rasters. A DNN or stack is not justified at this checkpoint: stacking
weak learners cannot correct the observed external-to-community level shift
without using the six local outcomes, and a more flexible DNN would add
overfitting risk rather than solve the transfer problem.

Route L remains the authoritative local mapping method. The external models and
their predictions are retained as a transparent negative result and as a guide
to data compatibility and future sampling.

## Reproducible files

- Preparation: `scripts/56_prepare_external_carbon_model.R`
- Model run: `scripts/57_fit_external_carbon_models.R`
- Functions: `R/covariates/external_carbon_model.R`
- Results: `outputs/covariate_analysis/part5_external_model/`

The result directory includes fold assignments and checks, every outer-fold
prediction, tuning results, per-community-core predictions, fitted model
objects, feature importance and observed-versus-predicted figures.

## Next step

Proceed to the Part 7 decision/reporting checkpoint rather than spatial
deployment. Part 6 flow-aware transfer remains deferred because raw TWI failed
the hydrology admission gate and the Part 5 residuals do not provide a
defensible basis for anisotropic mapping. A future external-model retry would
require a redesigned, locally comparable response set or substantially more
local complete profiles, not merely a more complex algorithm.
