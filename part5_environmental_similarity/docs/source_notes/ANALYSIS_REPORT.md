# Fort Severn hierarchical-similarity point test

## Outcome

The end-to-end point workflow runs correctly and the hierarchical similarity model improves on simple donor-source baselines. It is not yet accurate or well-supported enough to justify mapping the larger AOI.

The clearest issue is positive bias at 15–30 cm. FS-04 also has almost no close local analogue among the current Riley and CanPEAT donors and must be treated as an extrapolation case.

## Input audit

The input folder contains 12 CSVs: 11 donor shards and one community-validation shard.

| Check | Result |
|---|---:|
| Total points | 1,572 |
| External donor points | 1,564 |
| Community validation points | 8 |
| CanPEAT donors | 1,217 |
| Riley donors | 347 |
| Duplicate `core_uid` | 0 |
| Identical file schemas/order | Yes |
| Carbon/stock/coverage fields exported from Earth Engine | 0 |
| Covariate/QA bands | 47 |

Thirteen exported donor points—three CanPEAT and ten Riley—have no complete fixed-depth carbon interval in the local stock table. They are correctly excluded, leaving 1,551 usable donor cores.

The model uses 31 local, seven regional, and five broad continuous covariates. Four QA observation-count bands are retained for inspection but do not enter similarity distance.

## What the model does

1. Fit the median and MAD of every covariate using donor rows only. Use donor standard deviation only when MAD is zero.
2. Standardize donor and target values with those donor-only parameters.
3. Calculate separate local, regional, and broad similarities:

   \[
   S_s = \exp\left(-\frac{1}{2}\operatorname{mean}(z_{target}-z_{donor})^2\right)
   \]

4. Combine the three similarities using depth-specific exponents:

   \[
   K_d = S_L^{\alpha_{Ld}}S_R^{\alpha_{Rd}}S_B^{\alpha_{Bd}}
   \]

5. At each depth, use only donors with a complete carbon-stock interval.
6. Calculate one weighted prediction within CanPEAT and one within Riley.
7. Average the two source predictions equally so CanPEAT does not dominate merely because it has more cores.
8. Report prediction spread, effective donor counts, maximum similarity, block support, and top donors.

No community carbon measurements enter scaling, similarity, donor selection, or fitting. They are joined afterward for validation.

## Measured shallow validation

### All measurable community cores

| Interval | n | Mean observed | Mean predicted | Bias | MAE | RMSE | Correlation |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0–15 cm | 7 | 4.78 | 5.98 | +1.20 | 2.14 | 2.24 | 0.55 |
| 15–30 cm | 6 | 3.52 | 6.53 | +3.01 | 3.01 | 3.40 | 0.70 |
| 0–30 cm | 6 | 8.62 | 12.59 | +3.97 | 4.34 | 4.81 | 0.61 |

### Fort Severn cores only

| Interval | n | Mean observed | Mean predicted | Bias | MAE | RMSE | Correlation |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0–15 cm | 5 | 5.40 | 5.99 | +0.59 | 1.90 | 2.02 | 0.65 |
| 15–30 cm | 5 | 2.94 | 6.46 | +3.51 | 3.51 | 3.72 | 0.67 |
| 0–30 cm | 5 | 8.35 | 12.45 | +4.10 | 4.55 | 5.06 | 0.59 |

These correlations are descriptive only; five to seven validation observations are far too few for a stable performance claim.

![Observed versus predicted](community_validation_observed_vs_predicted.png)

## Does similarity help?

Yes, modestly. Compared with a constant baseline formed by averaging the CanPEAT and Riley donor means:

| Group | Interval | Similarity RMSE | Mean-baseline RMSE | RMSE improvement |
|---|---|---:|---:|---:|
| All community | 0–15 cm | 2.24 | 2.94 | 23.8% |
| FS only | 0–15 cm | 2.02 | 2.62 | 22.8% |
| All community | 15–30 cm | 3.40 | 4.17 | 18.4% |
| FS only | 15–30 cm | 3.72 | 4.55 | 18.3% |

The model is therefore learning useful environmental structure, but the donor response distribution remains too carbon-rich for these community cores—especially at 15–30 cm.

## Environmental support

Thirty-six of the 40 target-depth predictions pass the current maximum-similarity warning threshold. Four are low-support, all for FS-04 from 0–100 cm. Its maximum similarities are approximately 0.001, 0.003, 0.009, and 0.038 for the first four intervals. Its 100–200 cm value reaches only 0.17.

![Support heatmap](community_prediction_support_heatmap.png)

FS-04's local mismatch is driven most strongly by:

- extremely high Sentinel-1 VV/VH temporal-ratio variability;
- high Dynamic World water and flooded-vegetation probabilities;
- locally depressed elevation;
- unusually weak VH backscatter and low summer EVI.

This could represent a genuinely distinct wet microsite, a shoreline/water-contaminated pixel, a coordinate-placement issue, or an unstable radar-ratio metric. Inspect the point and its 30 m neighbourhood before changing model weights.

## Deep intervals

The model produced 30–50, 50–100, and 100–200 cm estimates, but the community data do not measure these intervals completely. They are predictions, not validations.

CanPEAT versus Riley disagreement increases with depth. Mean source predictions are about 49.6 versus 60.8 kg C/m² at 100–200 cm, with individual source gaps around 10.6–12.7 kg C/m². Deep outputs should therefore be labelled external-reference estimates until independent deep community cores are available.

## Decision

Continue developing the approach, but do not map the full AOI yet.

The next tests should be:

1. Inspect FS-04 against imagery, water masks, its coordinate, and the surrounding 3 × 3 pixels.
2. Test removal or robust capping of `local_s1_vv_over_vh_linear_temporal_sd`; ratios can become unstable near the VH noise floor.
3. Balance correlated local sensor families so several near-duplicate Dynamic World or Sentinel-1 bands do not collectively dominate the local block.
4. Test donor conditioning by peatland class/ecosystem and hydrologic setting.
5. Tune bandwidths and depth exponents only with held-out study/spatial folds, not these eight targets alone.
6. Add authoritative river/coast connectivity and surficial geology one family at a time, retaining the same validation gate.

## Reproduce the run

From the project root:

```bash
Rscript scripts/85_run_hierarchical_similarity_fs_point_test.R \
  --overwrite=true
```

This command audits inputs, runs the model self-test, fits the real point model, evaluates measured shallow intervals, regenerates the figures, and compares the result with simple donor baselines.
