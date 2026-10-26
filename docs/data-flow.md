# From spreadsheet to maps

Run code from the downloaded course R project, after `source("00_setup.R")`. Every workshop figure has copyable code directly below it. The code names the script that creates it and opens the exact saved image.

| Workshop step | Input | Run | Output used next |
| --- | --- | --- | --- |
| Part 1, Step 1 | [Workbook](../data/input/community_cores.xlsm): section measurements, core log, plot/site log | Import code in the [workshop](workshop.md#step-1-bring-field-and-lab-data) | Checked layers/core metadata and QC in memory. Step 2 rereads the same workbook. |
| Part 1, Step 2 and profiles | Same workbook; settings in [config.R](../config.R) | `source("scripts/01_field.R")` | [core_stocks.csv](../outputs/part1/core_stocks.csv), [layer_stocks.csv](../outputs/part1/layer_stocks.csv), [statistics](../outputs/part1/descriptive_statistics.csv), [QC](../outputs/part1/qc_flags.csv), three figures. |
| Part 1, Step 3 | Step 2 core/layer CSVs + [reference summaries](../data/reference/) | `source("scripts/02_compare.R")` | [Comparison figures/tables](../outputs/part1/comparison/). This is a context branch; it does not supply the spatial calculations. |
| Part 2A | Step 2 core CSV (stocks, ecosystem labels, coordinates) + [AOI GeoJSON](../data/aoi/study_area.geojson) | `source("scripts/03_area.R")` | [area_analysis.rds](../outputs/part2/area/area_analysis.rds): grid, area weights, AOI, core points/stocks. Also maps and totals. |
| Part 2B | Part 2A checkpoint | `source("scripts/04_spatial.R")` | [spatial_analysis.rds](../outputs/part2/spatial/spatial_analysis.rds): IDW values and a separate landscape simulation. Also maps and metrics. |
| Part 2C | Simulation, model map and test outputs saved by B | Read saved outputs; no second model run | Inspect simulation performance. This regression branch does not feed the Bayesian update. |
| Part 3 | **B's IDW** values in the spatial checkpoint + invented prior and assumed SDs | `source("scripts/05_bayes.R")` | [Updated summary](../outputs/part3/bayesian_summary.csv), [cell values](../outputs/part3/bayesian_grid.csv), map images and GeoTIFFs. |

The comparison branch adds reference context. The landscape branch demonstrates a model using additional simulated observations. The continuous workbook-to-update path is:

`workbook → core stocks → area checkpoint → IDW → Bayesian update`

The AOI boundary enters at Part 2A. Invented ecosystem areas, simulation inputs and Bayesian assumptions are identified in the [workshop](workshop.md) and each spatial folder's `method_metadata.csv`.

Each snippet requires its earlier output files. It does not automatically repair missing steps. For a complete rebuild after editing the workbook, use:

```r
source("run_all.R")
```

`R/view.R` only opens an existing saved figure in your browser/image viewer. It does not rerun an analysis. Run the creating script first to refresh a figure. The full course saves files without opening every image.

The source workbook is preserved. Step 2 records the confirmed example-only density correction and other QC checks alongside the recalculated outputs.
