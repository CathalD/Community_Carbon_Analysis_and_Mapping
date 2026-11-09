# Hierarchical similarity: direct Colab run card

## What this run produces

```text
External donor carbon table (local only)
             │ join later by core_uid
             ▼
Donor coordinates ─┐
                   ├─ Earth Engine: identical 47 covariates ─┐
Community points ──┘                                         │
                                                             ▼
                                      Colab /content → checked ZIP → browser
                                                             │
                                                             ▼
                                      R similarity model → FS validation
```

Carbon measurements, completeness flags, and response fields never enter Earth Engine.

## 1. Files to upload to Colab

- `scripts/gee/hierarchical_similarity_donor_sampling.py`
- `data/derived/external_carbon_profiles/external_carbon_cores_canada_for_gee.csv`
- `data/derived/community_cores_for_gee.csv`

## 2. First run: inspect only

Run the Python script unchanged. Confirm:

- exactly 47 bands pass server validation;
- the selected donor count and Riley/CanPEAT breakdown are plausible;
- the printed shard plan sums to the selected donor count;
- no export or download starts.

## 3. Second run: download point covariates

Use these settings in section 9:

```python
EXPORT_DESTINATION = "DIRECT"
RUN_POINT_EXPORTS = True
INCLUDE_COMMUNITY_VALIDATION_POINTS = True
RUN_DIRECT_RASTER_TILES = False
DOWNLOAD_DIRECT_ZIP_TO_BROWSER = True
```

Run the script again. Keep the browser tab open until this file is downloaded:

```text
HierarchicalSimilarity_DirectDownloads.zip
```

The ZIP contains point shards, the 47-band manifest, export plan, download manifest, hashes, and run metadata.

## 4. Optional raw AOI raster

Only after the point validation works, change:

```python
RUN_POINT_EXPORTS = False
RUN_DIRECT_RASTER_TILES = True
```

The full Fort Severn AOI is too large for one direct Earth Engine image request. The script therefore downloads aligned 4.5 km GeoTIFF tiles. Each tile is 150 × 150 pixels with 47 float32 bands, approximately 4.2 MB before compression.

For a single large GeoTIFF, keep using an Earth Engine batch export to Drive or Cloud Storage instead.

## 5. Unzip in this project

Place the ZIP contents here:

```text
data/raw/gee_imports/HierarchicalSimilarity_FS_AOI/
```

The model reads:

```text
data/raw/gee_imports/HierarchicalSimilarity_FS_AOI/point_shards/
```

Do not re-save the CSVs in Excel.

## 6. Run the model

From the project root:

```bash
Rscript scripts/83_run_hierarchical_similarity_micro_pilot.R --self-test

Rscript scripts/83_run_hierarchical_similarity_micro_pilot.R \
  --points=data/raw/gee_imports/HierarchicalSimilarity_FS_AOI/point_shards \
  --stocks=data/derived/external_carbon_profiles/external_carbon_cores_all.csv \
  --weights=config/hierarchical_similarity_depth_exponents.csv \
  --overwrite=true
```

## 7. Score the community and FS results

```bash
Rscript scripts/84_evaluate_hierarchical_similarity_community.R \
  --overwrite=true
```

Review these first:

- `community_validation_metrics.csv`: bias, MAE, RMSE, correlation, and mean support;
- `community_validation_rows.csv`: core-by-core residuals at measured intervals;
- `community_validation_0_30cm.csv`: complete shallow-profile comparison;
- `source_predictions.csv`: CanPEAT versus Riley disagreement;
- `top_donors.csv`: which profiles control each result;
- `block_support.csv`: local, regional, and broad similarity evidence.

## Decision gate

Do not map the larger AOI yet if FS error is poor, source disagreement is large, or maximum similarity/effective donor support is weak. First revise one covariate family or one prespecified depth-weight set, then repeat the same held-out test.
