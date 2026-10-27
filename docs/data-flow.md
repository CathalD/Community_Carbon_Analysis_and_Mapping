# From field sheet to figures

Start with [the five-step guide](start_here.md). Commands run from the downloaded R project. `source()` runs a saved script; the script creates the files named below.

## Main route

| Step | Input file | Run | Result | What happens next |
| --- | --- | --- | --- | --- |
| Try the example | Supplied workbook and boundary | `source("run_example.R")` | Five figures and tables in `outputs/example/` | Inspect the example, then prepare your two files. |
| Prepare your field sheet | `data/input/field_data_template.csv` | Replace all example rows; save as `data/input/my_field_data.csv` | One row for each measured slice | The next run checks depths, numbers and core details. |
| Add your boundary | Your polygon | Save as `data/aoi/my_aoi.geojson` | A boundary containing your sample locations | The next run uses it to limit the maps. |
| Run your data | Both files above | `source("run_my_data.R")` | Checked layers, stocks and three figures in `outputs/my_data/part1/` | Complete 0–30 cm stocks and coordinates feed the maps. |
| Read your maps | Checked stocks, coordinates and boundary | Created by the same run | Point map, distance-weighted map and GeoTIFF in `outputs/my_data/maps/` | Review the maps beside the samples and data checks. |

**The main chain is:**

`my_field_data.csv → checked slices → core_stocks.csv`

`core_stocks.csv + my_aoi.geojson → sample_map.png and idw_stock_map.png`

`my_aoi.geojson` enters at the mapping step. It sets the map boundary; it does not supply carbon values. Replacing the boundary alone does not move sample locations.

## Copy code instead of running the wrapper

This does the same main analysis using your files:

```r
source("00_setup.R")
source("R/course.R")
run_course(
  workbook_path = "data/input/my_field_data.csv",
  aoi_path = "data/aoi/my_aoi.geojson"
)
```

The default output folder is `outputs/my_data/`. The name `workbook_path` also accepts the simple CSV field sheet. You do not need to change calculation functions.

Every main-route image in the [workshop](workshop.md) has a short block below it that runs the analysis and opens that image. Running the analysis again refreshes results; opening an image alone does not.

## Optional full illustrative course

The extensions retain their own example inputs in `config.R` and outputs in `outputs/part1/`, `outputs/part2/` and `outputs/part3/`. They do not use or overwrite `outputs/my_data/`.

```r
source("run_all.R")
```

| Optional step | Input | Run one step | Result |
| --- | --- | --- | --- |
| Example field analysis | Supplied workbook via `config.R` | `source("scripts/01_field.R")` | Stocks and profiles in `outputs/part1/`. |
| Reference comparisons | Example core/layer stocks + supplied reference summaries | `source("scripts/02_compare.R")` | Context figures in `outputs/part1/comparison/`. This is a separate branch. |
| Ecosystem totals | Example stocks + `data/aoi/study_area.geojson` | `source("scripts/03_area.R")` | Illustrative zones, totals and saved area results in `outputs/part2/area/`. |
| Spatial examples | Saved area results | `source("scripts/04_spatial.R")` | IDW map and separate landscape simulation in `outputs/part2/spatial/`. |
| Bayesian update | IDW results + invented prior and assumed uncertainty | `source("scripts/05_bayes.R")` | Updated example maps in `outputs/part3/`. |

These extensions demonstrate additional ideas. The invented ecosystem zones, landscape observations and prior map are not real inputs for your own AOI. See the [full workshop](workshop.md#optional-full-illustrative-course) and [source details](../SOURCES.md).
