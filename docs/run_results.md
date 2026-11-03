# Full four-method rerun

The current comparison has **four local-method branches and four prior updates for each analysis**: Sothe S30, Li organic depth and Li organic carbon. This is 12 updated raster products, excluding the three original priors. `outputs/verification/completed_updates.csv` lists each product. Current figures are under `outputs/results/full_suite/`. Earlier two-contact results and their figures remain historical.

All eight FS cores are retained. Sothe uses the unchanged target-or-sampled-bottom stock function. Li uses all eight contact estimates authorized by the user, with two confirmed and six layer-classified contacts separately labelled. Four layer-classified responses are zero organic depth/stock. No Route L mapping is run.

## Mean mapped estimates

| Analysis / units | Prior | Spatial core update | Area update | IDW update | Regression update |
|---|---:|---:|---:|---:|---:|
| Sothe S30, kg C/m² | 31.39 | 29.98 | 13.18 | 13.99 | 13.74 |
| Li depth, cm | 125.30 | 121.19 | 22.12 | 10.66 | 26.15 |
| Li organic carbon, kg C/m² | 64.14 | 62.20 | 3.12 | 1.28 | 10.70 |

Each row uses a common valid prior footprint across its four methods: Sothe 3,873 cells; Li 3,831 cells. Regression cells without physical local support retain the prior. The differences describe alternative revisions to estimates, not measured carbon loss. Pointwise map means are not uncertainty estimates for whole-area inventories.

Sothe retains the workshop log-Gaussian method and the S100-derived S30 SD sensitivity proxy. Script 3 selects Mean_only; its fourth regression update is therefore intercept-only, not a forced covariate model.

## Li scale adaptation and support

Li uses original-scale Gaussian conditioning/fusion so authorized zero responses can enter without an arbitrary offset. This is explicitly different from the workshop's positive lognormal stock assumption. Gaussian ensemble-spread uncertainty is a sensitivity, not calibrated predictive error.

The spatial updates produce **33 negative depth means and 14 negative carbon means**. These are flagged and retained for diagnosis, not clipped or interpreted as physical negative stock. Local regression also extrapolates below zero; those estimates are not fused as evidence. The prior is retained in **359 depth cells and 476 carbon cells**. These limits should accompany the maps.

Li source NoData at FS-04 and FS-02 requires paired conditioning priors 18.96 m and 21.85 m away, respectively, with a 500 m maximum. All extractions are audited. Valid zero-valued source depth cells are preserved.

Both Li regressions select Shrubs under the workshop AICc rule. The 119 bands are audited; models use a prespecified small subset. Fixed-selection local errors are not the procedure-validation errors below.

## Outer whole-core validation of the Li updates

IDW tuning, regression selection and working-error estimation are refit inside every outer training fold. Values are RMSE in the corresponding row's units.

| Variable | Prior | Spatial update | Area update | IDW update | Regression update |
|---|---:|---:|---:|---:|---:|
| Organic depth, cm | 102.71 | 76.38 | 23.54 | 19.25 | 70.67 |
| Organic carbon, kg C/m² | 56.09 | 39.09 | 4.05 | 4.19 | 34.08 |

These are eight nearby cores, six with uncertain contact classifications. They do not establish independent-site performance or calibrated uncertainty. Full fold predictions and model choices are saved beside each update. In this run, the covariate update does not outperform the simpler area/IDW updates in outer validation.

## Carbon-depth curves

Four positive organic profiles provide the cumulative-carbon slope estimate **0.2861 kg C/m² per cm**, with three residual degrees of freedom and positive-core holdout RMSE **2.2260 kg C/m²**. The four zero-depth profiles remain displayed at the origin and do not inflate curve precision. The deepest measured organic contact is 30.73 cm.

Every updated depth method has a separate `curve_derived_carbon.tif`, evaluated only within measured depth support. This supplies a depth-to-carbon diagnostic for each branch. The direct Li carbon prior updates remain separate marginal analyses; there is no invented deep tail, independent fusion of Li depth/carbon priors, or complete joint depth/curve uncertainty propagation.

More confirmed contacts, deeper complete peat profiles, and field/lab error measurements would improve this analysis. Temporal accumulation rates require chronology. Numerical verification checks confirm the implemented calculations and branch completeness; they do not establish ecological validity of the contact classifications or Gaussian-scale assumptions.
