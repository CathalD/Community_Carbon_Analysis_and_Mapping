# Community Carbon Stock Analysis and Mapping

Turn a soil field sheet into carbon-stock tables, depth profiles and a simple map.

**Start with the supplied example. Then use your own measurements and area boundary.** The same steps work with a local study area in Canada when you supply suitable data.

[Start here: every file and command](docs/start_here.md) · [Workshop: figures and code](docs/workshop.md) · [Explainer Google Doc](https://docs.google.com/document/d/1qflFw6CH69feNRr87Oi1SbHl0sY5bcXb-TuFI6KUtaU/edit)

## Download and run the example

In R or RStudio, paste:

```r
install.packages("usethis") # Once only.
usethis::use_course(
  "https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip"
)
```

Open `Community_Carbon_Analysis_and_Mapping.Rproj` in the downloaded folder. Paste this into the **Console**:

```r
source("run_example.R")
```

`source()` runs the saved R script. The script installs missing packages, reads the example workbook and makes five figures. Find them in `outputs/example/`. The first run needs internet access to install packages.

## Use your data

1. Copy `data/input/field_data_template.csv` to `data/input/my_field_data.csv`. Open it in a spreadsheet program, replace **all** example rows, keep the eight column headers, and save as CSV.
2. Save your study boundary as `data/aoi/my_aoi.geojson`.
3. Run:

```r
source("run_my_data.R")
```

Your tables and five figures go to `outputs/my_data/`. The [start guide](docs/start_here.md) explains the columns, boundary and output files. The [workshop](docs/workshop.md) has code directly below each figure.

Reference comparisons, ecosystem totals, landscape models and Bayesian updating are [optional extensions](docs/workshop.md#optional-full-illustrative-course). Course code is [MIT licensed](LICENSE); see [source details](SOURCES.md) for third-party data.
