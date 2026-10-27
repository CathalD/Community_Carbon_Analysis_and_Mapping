# Start here

This guide prepares the files and software for the complete workshop. The [workshop](workshop.md) explains each analysis, shows its example outputs and places runnable code below the figures. Work through the sections in order.

**Workshop sequence:** field measurements → carbon stocks and descriptive statistics → reference comparisons → ecosystem area estimates → spatial estimates and model assessment → Bayesian updating.

**You need:** [R](https://cran.r-project.org/) and [RStudio Desktop](https://docs.posit.co/ide/user/download.html). Install R first, then RStudio. The supplied example needs no file changes.

## 1 Download and open the course

Paste this into the RStudio **Console**:

```r
install.packages("usethis") # Once only.
usethis::use_course(
  "https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip"
)
```

Open `Community_Carbon_Analysis_and_Mapping.Rproj` in the downloaded folder. Keep the folder structure unchanged. Opening the project makes the course folder the starting point for every file path in the code.

## 2 Run the complete supplied example

```r
source("run_all.R")
```

`source()` runs a saved R script. This command runs all five analysis scripts in order. Setup installs missing packages: `readxl` reads spreadsheets, `ggplot2` draws figures, `sf` handles boundaries and distances, and `terra` writes map files. Calculations, regression and Bayesian updating also use base R and `stats`, included with R. The first run needs internet access. Wait for installation to finish.

The example uses four illustrative cores in western James Bay. Later stages include hypothetical ecosystem zones, simulated landscape data and an invented prior map. These inputs let you follow the complete analysis; they do not establish measured stocks or map accuracy for a real community.

## 3 Work through the analyses in order

To run one stage at a time, first run setup:

```r
source("00_setup.R")
```

Then use the commands in this table, in order. They reproduce the same results as `run_all.R`.

| Workshop section | What you analyse | Command | Results folder |
| --- | --- | --- | --- |
| Part 1, Steps 1–2 | Field/lab data checks; layer and core stocks; depth profiles; means, standard deviations and coefficients of variation. | `source("scripts/01_field.R")` | `outputs/part1/` |
| Part 1, Step 3 | Stock, carbon concentration and density relative to reference collections, studies, ecosystems and regions. | `source("scripts/02_compare.R")` | `outputs/part1/comparison/` |
| Part 2A | Ecosystem areas, mean stock maps and area totals. | `source("scripts/03_area.R")` | `outputs/part2/area/` |
| Parts 2B–C | Inverse-distance weighting (IDW), distance to samples, a landscape regression and held-out prediction checks. | `source("scripts/04_spatial.R")` | `outputs/part2/spatial/` |
| Part 3 | Prior, new and updated carbon estimates; assumed uncertainty. | `source("scripts/05_bayes.R")` | `outputs/part3/` |

The reference comparison uses the field outputs. The area results feed the spatial stage. Bayesian updating uses the IDW estimate from Part 2B; the landscape model in Part 2C is a separate demonstration.

In the modelling stage, **held-out data** are observations kept out of model fitting so that you can compare predictions with values the fitted model has not used.

Open figures from their results folder, or run the figure-specific code in the [workshop](workshop.md). For example:

```r
source("R/view.R")
view_figure("outputs/part1/depth_profiles.png")
```

After changing an input, rerun its analysis and the stages that use its results. See [how the files connect](data-flow.md).

## 4 Prepare your own field measurements

After completing the supplied analysis, practise with your own measurements. Copy [field_data_template.csv](../data/input/field_data_template.csv) and save it as **`data/input/my_field_data.csv`**.

Open the copy in a spreadsheet program. Replace **all** example rows. Keep the headers unchanged. Use one row per sampled slice.

| Column | What to enter |
| --- | --- |
| `core_id` | A name such as `CORE-01`; repeat it for every slice from that core. |
| `ecosystem` | A label such as `Bog` or `Forest`; repeat it for that core. |
| `latitude` | WGS84 decimal degrees, such as `52.6`; repeat the core's location. |
| `longitude` | WGS84 decimal degrees, such as `-84.3`; west is negative. |
| `top_cm` | Top of the slice, in centimetres below the ground surface. |
| `bottom_cm` | Bottom of the slice, using the same surface. |
| `bulk_density_g_cm3` | Dry bulk density in g/cm³, such as `0.09`. |
| `organic_carbon_pct` | Organic carbon as a percent: enter `45` for 45%. |

Start depth at **0 at the ground surface**, including the organic surface layer. Slices must not overlap. A 0–30 cm stock needs complete coverage from 0 to 30 cm; a gap or a shallower core gives a missing value for that depth. Keep numbers numeric: enter `10`, not `10 cm`.

Save as a comma-separated **CSV** file. Alternatively, put this table on the first sheet of an `.xlsx` file and edit the `my_field_data` path at the top of `run_my_data.R`. See the [column guide](../data/input/field_columns.md) for an example table.

## 5 Add your area of interest

An **area of interest (AOI)** is the boundary of the place you want to map. Save a polygon boundary as **`data/aoi/my_aoi.geojson`**. GeoJSON stores the outline and its geographic position. Export an existing boundary with a known coordinate system; WGS84 / EPSG:4326 uses longitude and latitude in degrees.

To draw a boundary, open [geojson.io](https://geojson.io/) in your browser. Zoom to your study area, choose the polygon tool, click the corners and close the outline. Download/export GeoJSON, rename it `my_aoi.geojson`, and move it to `data/aoi/`. See the [tool documentation](https://github.com/mapbox/geojson.io).

Use a boundary containing your samples. Changing the boundary does not move the example cores: replace their coordinates too. The IDW code requires at least **two different locations inside the boundary**, each with a complete 0–30 cm stock. This is a software minimum; two cores do not establish representative sampling or map accuracy.

## 6 Practise field calculations and IDW with your data

```r
source("run_my_data.R")
```

With the filenames above, no code changes are needed. This practice script runs field calculations and IDW, saving to `outputs/my_data/`. It preserves your input files. `run_example.R` runs the same practice sequence with the supplied data and saves to `outputs/example/`; use `run_all.R` for the complete workshop.

| File inside `outputs/my_data/` | What to inspect |
| --- | --- |
| `part1/qc_flags.csv` | Data checks to review first. |
| `part1/core_stocks.csv` and `part1/layer_stocks.csv` | Core and slice calculations. |
| `part1/descriptive_statistics.csv` | Group means and spread. |
| `part1/common_depth_stock.png` | Stocks compared over 0–30 cm. |
| `part1/depth_profiles.png` and `part1/cumulative_stock.png` | Carbon and density down each core. |
| `maps/sample_map.png` | Boundary and usable sample locations. |
| `maps/idw_stock_map.png` and `maps/idw_stock_0_30.tif` | Estimated stocks between samples. |
| `maps/excluded_cores.csv` and `maps/mapping_metadata.csv` | Exclusions and mapping settings. |

To adapt the later analyses to a real area, you also need representative ecosystem areas, suitable landscape predictors and observations for model assessment, and a prior estimate with justified uncertainty. The teaching scripts generate these inputs for the supplied example. Substituting an AOI alone does not supply them. The [data-flow guide](data-flow.md) identifies these requirements.

## If something stops

- **File not found:** open the course `.Rproj`; check the exact filenames and folders.
- **Missing column or invalid number:** compare your field sheet with the headers and units above.
- **No mapped samples:** check longitude signs, the boundary and complete 0–30 cm coverage.
- **Package installation failed:** read the installation message. Some systems need extra system libraries for `sf` or `terra`; ask your instructor or computer support for help.

## Continue with real local analyses

The complete teaching example remains `source("run_all.R")`, including ordinary kriging on separate simulated observations. To add reference comparisons and optional real ecosystem, regression and prior inputs to your own-data exercise, use `source("run_local_workshop.R")`. Follow [the exact local interfaces](local-inputs.md). Missing local inputs do not trigger simulated substitutes.
