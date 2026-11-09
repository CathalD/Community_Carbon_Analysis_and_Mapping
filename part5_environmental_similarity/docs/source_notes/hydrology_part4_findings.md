# Part 4 hydrology findings

**Completed:** 25 September 2026  
**Status:** processing valid; raw TWI not admitted to modelling

## Processing verification

The downloaded Earth Engine cloud-optimized GeoTIFFs were rewritten as plain
Float32 GeoTIFFs with explicit NoData before WhiteboxTools processing. The
rewritten CDEM and Copernicus inputs match their downloads in geometry,
valid-cell count, minimum, mean and maximum. Both candidates completed
least-cost breaching, depression/flat resolution, slope, D-infinity flow
direction, specific catchment area and strict TWI without a processing error.

Both candidates cover all 4,802,814 AOI cells. Strict TWI is missing at only 14
CDEM cells and 10 Copernicus cells. Earlier summaries reporting approximately
50% missingness counted the area outside the irregular AOI but inside its
rectangular raster extent; the summary code and files have been corrected.

## DEM conditioning

| Check | CDEM | Copernicus |
|---|---:|---:|
| Absolute conditioning change >0.1 m | 6.59% | 6.29% |
| Absolute conditioning change >1 m | 3.05% | 0.69% |
| Absolute conditioning change >5 m | 0.014% | 0.008% |
| Slope <0.1 degrees | 9.27% | 14.60% |

Most changes are small, but CDEM required materially more adjustment above
1 m. Its integer elevations also create extensive flats. Copernicus preserves
continuous elevations but is a surface model and includes hydro-edited water.

## Candidate agreement

The conditioned elevation surfaces agree closely across the AOI (Pearson
correlation 0.998; Copernicus averages 2.51 m above CDEM). Hydrologic derivatives
are substantially less stable:

| Layer | AOI Pearson correlation |
|---|---:|
| Slope | 0.608 |
| Specific catchment area | 0.157 |
| Strict TWI | 0.679 |

At the eight community cores, the correlation is 0.172 for slope, -0.213 for
log10 specific catchment area and -0.362 for TWI. The median absolute TWI
difference is 1.80 and the maximum is 6.96. FS-01 and FS-03 differ by almost
seven TWI units, with particularly large disagreement in contributing area.

Strict TWI also has a strong near-zero-slope tail. TWI exceeds 30 over 5.47% of
the CDEM AOI and 3.73% of the Copernicus AOI. These values are mathematically
consistent with the strict formula but are highly sensitive to small slope and
conditioning differences in this low-relief landscape.

## Copernicus quality masks

Copernicus marks 93.04% of AOI cells as unedited and 93.04% as not filled.
Water covers 6.80% of the AOI: 3.64% ocean, 1.90% lake and 1.26% river. The HEM
is present for 93.04% of cells, with median 0.373 m and 95th percentile 0.581 m.
All eight community cores are unedited, unfilled, non-water pixels. Their HEM
values range from 0.325 to 1.739 m.

## Decision

The Whitebox workflow itself passes technical QA, but raw strict TWI does not
pass the scientific admission gate because the routed quantities are too
DEM-sensitive at the community cores and the extreme tail is strongly driven
by near-zero slopes. Do not add either raw TWI raster to the carbon model yet.

Retain the existing slope and relative-elevation covariates. If hydrology is
pursued further, the next bounded test should validate mapped flow lines
against independent hydrography and compare a transparently regularized or
rank-based wetness representation across both DEMs. Any such layer remains a
sensitivity candidate until it is stable at the cores and improves held-out
carbon or organic-depth performance.

## Audit outputs

- `outputs/covariate_analysis/part4_hydrology/hydrology_candidate_qa_metrics.csv`
- `outputs/covariate_analysis/part4_hydrology/community_core_hydrology_comparison.csv`
- `outputs/covariate_analysis/part4_hydrology/candidate_comparison.csv`
- Candidate-specific summaries, per-core tables and rasters under `cdem/` and
  `copernicus/`

