# Carbon Analysis and Mapping Workshop

Open `Carbon_Analysis_and_Mapping_Workshop.Rproj` in RStudio. Run the numbered scripts from this project folder. All paths are relative to this folder, so the project can be moved without editing an absolute computer path.

## Add your inputs

- **Field workbook:** `data/field/field_cores.xlsm`. The supplied workbook has six simulated cores. Replace it with a workbook using the same sheet names and column headings, or change `workbook_path` near the top of `script_1.R`.
- **Study boundary:** `data/boundary/aoi.geojson`. Use a correctly identified coordinate system and field coordinates belonging to that study area.
- **Reference observations:** `data/reference/`. The existing five source collections, summaries, ecosystem stock CSV, source links and provenance are included. Script 1 reads these local files without downloading data.
- **Environmental layers:** covariate raster folders are excluded from this branch. Before running Script 3, supply your local rasters as `data/covariates/topography.tif` and `data/covariates/water.tif`. The covariate export code remains in `tools/`.
- **Prior maps:** `data/priors/`. The supplied McMaster stock and uncertainty rasters are existing prepared 1 km files. They are not the original native-resolution releases.

## Run the workflow

| Script | Reads | Writes |
| --- | --- | --- |
| `script_1.R` | Field workbook and reference observations | `outputs/script1/`: mapping table, core summary, profile plots, reference comparisons and Part 3 summary |
| `script_2.R` | Script 1 mapping table and study boundary | `outputs/script2/`: area-based polygons, study boundary, IDW maps, tuning and totals |
| `script_3.R` | Script 1 mapping table, study boundary and covariate rasters | `outputs/script3/<stock variable>/`: extracted predictors, model fits, validation and maps |
| `script_4.R` | Script 1 table, Script 2 maps and supplied priors | `outputs/script4/`: comparison table and aligned prior rasters |
| `script_5.R` | Script 1 table, Script 2 maps and Script 4 prior rasters | `outputs/script5/`: Bayesian model products and validation |

Start with Script 1 to inspect the field measurements. Continue with Script 2 for area averages and IDW mapping. The later scripts keep their existing analysis code for the remaining workshop work.

`outputs/script1/mapping_data.csv` is the single field-stock table used by the mapping scripts. Do not manually edit a generated output to change an observation. Update the workbook and rerun Script 1, then rerun the affected later scripts.

Every script creates its own output directory. Results from a rerun replace that script's generated files. The current included results use the simulated example observations; they are not a measured community inventory.

## Stock depth groups

The agreed rules are unchanged:

- **0–30 cm:** cores reaching 30 cm, capped at 30 cm.
- **0–50 cm:** cores deeper than 30 cm, capped at 50 cm or their sampled bottom.
- **0–100 cm:** cores deeper than 50 cm, capped at 100 cm or their sampled bottom.
- **Total:** the measured sampled profile.

Sections crossing a reporting depth contribute their stock per centimetre multiplied by the centimetres within the window. No material below the sampled bottom is extrapolated. Partial target-depth coverage remains identified in the reference data. Sampled profile totals span different depths and do not establish mineral contact or the entire soil carbon store.

## Other folders

- `tools/export_covariates.ipynb`: the existing Earth Engine/Colab covariate workflow. Set its project, dates and export settings when using it.
- `tools/prepare_reference.py`: rebuild the four reference depth groups and their summary tables from the bundled prepared layers.
- `tools/export_ecosystem_stocks.py`: rebuild the combined ecosystem stock CSV after updating the reference groups.
- `helpers/map_helper.R`: optional context-map functions. These functions expect `data/context/north_america.gpkg` when called; that context file was not part of this transfer.
- `docs/`: workflow information, file-name mapping, checks and the Script 6 overview.
- `archive/transfer_original/`: an unchanged snapshot of the supplied transfer, including its earlier outputs, duplicate file names and paste guides. The active scripts do not read from the archive.

To refresh the reference data preparation from the project folder:

```text
python3 tools/prepare_reference.py
python3 tools/export_ecosystem_stocks.py
```

Third-party source terms and attribution remain in `data/reference/resources.csv` and the reference README. No GitHub repository was changed by this local reorganization.
