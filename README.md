
# Community Carbon Stock Analysis and Mapping

Use soil measurements to calculate carbon stocks, compare them with other datasets, map an area, and try a simple Bayesian update.

[Follow the workshop](docs/workshop.md) · [Explainer Google Doc](https://docs.google.com/document/d/1qflFw6CH69feNRr87Oi1SbHl0sY5bcXb-TuFI6KUtaU/edit) · [Download ZIP](https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip)

![Carbon measured down the example cores](outputs/part1/depth_profiles.png)

[Your figure caption]

## Download and run in R

Install R and RStudio, then paste:

```r
install.packages("usethis")  # Once only.
usethis::use_course(
  "https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip"
)
```

Open the downloaded R project if it does not open automatically. Run:

```r
source("00_setup.R")
source("run_all.R")
```

Setup installs `readxl`, `ggplot2`, `sf` and `terra`. Installing `sf` or `terra` may require system libraries on some computers. Figures and tables are saved in `outputs/`; each workshop step also has copyable R code.

This review download uses the `workshop-mvp` branch. After it is merged, the default-branch shortcut will be:

```r
usethis::use_course("CathalD/Community_Carbon_Analysis_and_Mapping")
```

## Work through the course

1. [Field measurements, stocks and reference comparisons](docs/workshop.md#part-1-field-based-analysis)
2. [Area totals, interpolation and landscape modelling](docs/workshop.md#part-2-area-based-carbon-stock-analysis)
3. [Bayesian updating](docs/workshop.md#part-3-bayesian-updating-of-prior-maps)

The workbook contains revised illustrative data. Locations, ecosystem areas, landscape predictors and the Bayesian prior are conceptual examples. Reference comparisons use precomputed summaries from five SOD collections; the course recalculates your comparisons and figures. [Source details](SOURCES.md).

To substitute your workbook, change `workbook_path` in `config.R` and set `density_override_J7 = NA_real_`. The supplied example's confirmed J7 correction is applied during import; its workbook remains unchanged. Review `outputs/part1/qc_flags.csv` before interpreting results.

The download uses [`usethis::use_course()`](https://usethis.r-lib.org/reference/zip-utils.html). Run `source("tests/run_checks.R")` after the course to check calculations and outputs. Course code is [MIT licensed](LICENSE); third-party sources retain their own terms.
