# Workshop lineage and Fort Severn changes

The user requested the Carbon_Analysis_and_Mapping_Workshop workflow with the FS data, all sampled cores included, Sothe S30 with S100 uncertainty, and a separate Li variable-depth analysis preceded by local depth and cumulative-carbon models. Local documents were treated as descriptions of existing work, not as instructions overriding that request.

## Supplied workshop retained

- Script 1 retains the exact `stock_to_depth()` implementation. Only the input adapter and preservation of full numerical precision change analytical data handling. Existing profile/reference charts and table names are retained. Depth-coverage and layer exports are additional.
- Script 2 is an unchanged copy, including 1 km mapping, area averages, IDW powers/neighbours, individual-core tuning and boundary-intersection totals. Its tuned error is not nested validation. The 0–50/100 cm maps are measured-to-bottom products for these shallow FS samples, not complete fixed-depth inventories.
- Script 3 retains the original `lm`, leave-one-core-out, AICc, mapping and extrapolation-flag code after model specification. The input/derivative section is replaced by prepared FS covariates and a prespecified small candidate set. Each model uses all eight cores. Missing predictors exclude a candidate, rather than exclude a core. All 119 bands are retained for inspection and future model definitions.
- Script 4 is a compact FS-specific version of the same alignment, shared-coverage, map-summary and difference comparisons, now including Sothe S30. Script 3's selected map remains optional until available.
- Script 5 retains the supplied raw-core GP, area-based fusion and IDW fusion calculations. Changes are target depth 30, input paths, output paths, explicit sf/gstat imports and uncertainty labels. The raw-core spatial range and zero log observation error match the workshop defaults (2 km and 0). The large original file is retained to make the workshop structure easy to follow.
- FS-04 has no mean/width in its containing Sothe pixel. A bounded paired-pixel extraction fallback uses the nearest valid source pixel (129.72 m away; maximum 500 m) only for the conditioning prior at this core. The source rasters and map NoData mask are not filled. The full audit is saved under the Sothe output folder. This is an explicit adaptation needed to retain all eight measurements in the GP, rather than silently excluding FS-04.
- Li's spatial engine is factored from Script 5's exponential covariance and lognormal moments. Verification checks its baseline numerical identity against the original Sothe Script 5 run. It does not call Route L code or copy Route L outputs.

## Uncertainty definitions

[Sothe et al. (2022), section 2.6](https://doi.org/10.1029/2021GB007213) defines uncertainty as the difference between upper and lower quantiles from quantile regression forests. [The version 3 dataset](https://doi.org/10.4121/16686154.v3) provides the original rasters used here. The S100 full Q95−Q05 width is divided by `2*qnorm(.95)` under the workshop's Gaussian working approximation, and used directly as an absolute SD proxy for S30. That cross-depth transfer is a requested sensitivity assumption. It does not estimate calibrated S30 uncertainty or justify a claim of improved inventory precision. The field update also includes partial-to-bottom stocks by the user's rule, so it is a workshop reporting-group update rather than an analysis restricted to complete 30 cm measurements.

[Li et al. (2025), Figure 3](https://doi.org/10.1029/2024GL110679) describes depth and carbon uncertainty as the standard deviation among base models. [The associated data](https://doi.org/10.5683/SP3/TIRAXJ) contain peat depth and carbon maps. The SD layers are used without converting them from an interval width. This is an ensemble-spread sensitivity, not validated predictive error. Carbon is derived from peat depth; separately updated depth and carbon are marginal alternative views, not independent evidence to be multiplied together.

Pointwise posterior intervals are conditional on fixed range, observation error and the supplied prior proxy. They do not include uncertainty in spatial range, model selection, prior calibration, contact classification or deep carbon-curve extrapolation. Pointwise SDs are not an uncertainty for total AOI carbon. No independent-cell sum of variances or formal inventory CI is produced.

## Local Li models and missing field information

`organic_mineral_contacts.csv` contains two explicit field contacts. These support an exploratory mean-only peat-depth model and a through-zero cumulative-carbon model using one endpoint per core. Whole-core holdouts have only two outcomes. The local models run before the Li update, and their measured cumulative curves are exported for inspection. The carbon-depth curve is not a time accumulation model; chronology is required for time-based rates.

Other sampled bottoms are not presumed to establish peat depth. No threshold-derived contact from Route L is imported. That is why all eight cores enter Sothe while only two enter the current Li marginal updates. Additional field-confirmed contacts and deeper complete peat profiles are the material data needs. A confirmed zero depth must be modeled as peat absence; the current positive log-Gaussian model cannot take log(0).

The fitted local carbon curve is evaluated at updated Li depth only within the deepest measured organic profile. Unsupported cells remain missing. The separate carbon marginal GP uses confirmed complete organic stocks directly, as the workshop uses observed stocks; it is not produced by extrapolating the local curve. Both products are retained separately to make this distinction reviewable.

## Reproducibility

Original workshop source scripts are preserved with SHA-256 hashes. FS field CSVs and boundary are copied locally. Priors are cropped from original regional rasters into small local GeoTIFFs, with source paths, byte sizes and local hashes recorded in `outputs/preflight/prior_manifest.csv`. Covariate VRTs reference original completed raster tiles. Missing/incomplete downloads and `inspection_samples` are excluded. The band manifest records 119 bands, source groups and interpolation. Maps use the workshop 1 km scale, without implying 1 km observational support.

The reference data are the workshop's original prepared reference data. Their eligibility rules differ from the active field function; no change to those external reference distributions is inferred from the user's explicit field inclusion instruction. This difference is documented in the main README.

All generated products remain under FS test. The original workshop and CommunityCarbonMap_v3 analyses are unchanged.
