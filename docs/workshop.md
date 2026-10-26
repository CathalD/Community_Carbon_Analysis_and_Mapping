# Community Carbon Stock Analysis and Mapping

This workshop shows what we analyse, how to read the results, and how to use them. It follows the [explainer Google Doc](https://docs.google.com/document/d/1qflFw6CH69feNRr87Oi1SbHl0sY5bcXb-TuFI6KUtaU/edit).

Run the example first, inspect its outputs, then substitute your data. Start with the [download and setup instructions](../README.md#download-and-run-in-r).

Work through these handoffs in order. The [input/output guide](data-flow.md) lists the files used at each step.

1. Part 1, Steps 1–2: workbook → checked layers → core stocks and profiles.
2. Part 1, Step 3: core stocks/layers + reference summaries → comparisons (a context branch).
3. Part 2A: core stocks + area boundary → ecosystem areas and totals.
4. Part 2B: area checkpoint → distance-weighted map and a saved spatial checkpoint.
5. Part 2C: inspect a separate landscape simulation saved in that checkpoint.
6. Part 3: use Part 2B's distance-weighted map + an invented prior → updated map.

Run setup once. Then paste the code below each figure from the course R project. Each block names its input files; run the earlier steps first. `view_figure()` opens the saved image in your computer's browser or image viewer. It does not recalculate it.

The revised workbook is illustrative. The 350 km² western James Bay area, core locations and ecosystem assignments are conceptual. Later landscape and Bayesian examples use additional simulated inputs, clearly separated from the four cores.

## Part 1 Field based analysis

### Step 1 Bring field and lab data

**Input:** [Workbook](../data/input/community_cores.xlsm): `3. Peat Data`, `2. Core Log`, and `1. Plot & Site Log`.

**Output:** Checked `core_input` in R, including QC. **Used next:** Step 2 recalculates stocks from the same workbook.

**Analyse:** Bring each core's depth intervals, dry bulk density and organic carbon into [one workbook](../data/input/community_cores.xlsm). `readxl` reads the inputs; macros are not needed.

**Read:** Check IDs, units and depth coverage. The code below shows QC in memory; Step 2 saves [qc_flags.csv](../outputs/part1/qc_flags.csv). The confirmed J7 density correction, 0.9 → 0.09 g/cm³, is applied during import; the source workbook stays unchanged.

**Apply:** Correct data-entry issues before using the results. For another workbook, change `config.R` and disable the example-only correction.

<details>
<summary>Copy R code</summary>

```r
source("config.R")
source("R/field.R")
core_input <- import_core_workbook(
  course_config$workbook_path,
  density_override_J7 = course_config$density_override_J7
)
core_input$qc
```

[Import and checking methods](../R/field.R)

</details>

### Step 2 Calculate stocks

**Input:** The [workbook](../data/input/community_cores.xlsm) checked in Step 1.

**Output:** [Core stocks](../outputs/part1/core_stocks.csv), [layer stocks](../outputs/part1/layer_stocks.csv), [statistics](../outputs/part1/descriptive_statistics.csv), [QC](../outputs/part1/qc_flags.csv), and the three figures below. **Used next:** Step 3 reads core/layer stocks; Part 2A reads core stocks and coordinates.

![Stocks over the same depth](../outputs/part1/common_depth_stock.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Re-read the workbook and rebuild Step 2 outputs.
source("scripts/01_field.R")
source("R/view.R")
view_figure("outputs/part1/common_depth_stock.png")
```

</details>

**Analyse:** Calculate carbon in each slice, then add slices over a shared depth:

`Stock (kg C/m²) = 0.1 × bulk density (g/cm³) × organic C (%) × thickness (cm)`

Enter 45% as `45`. Multiply kg C/m² by 10 for t C/ha.

**Read:** The revised 0–30 cm stocks are 14.82 and 35.20 kg C/m² for the forest cores, and 10.44 and 12.70 for the wetland cores. Deeper sampled totals cover different depths. Forest 0–100 cm values remain missing.

**Apply:** Compare equal depths. Missing soil below a core is not zero carbon.

[Stock calculations and figure code](../R/field.R)

### Step 2 continued Read carbon down the core

**Input:** The same [workbook](../data/input/community_cores.xlsm). Each block below reruns Step 2 and opens its named figure.

**Output:** [Layer stocks](../outputs/part1/layer_stocks.csv) and profile figures. **Used next:** Step 3 uses the layer table to match reference depths.

![Carbon and density profiles](../outputs/part1/depth_profiles.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Re-read the workbook and rebuild Step 2 outputs.
source("scripts/01_field.R")
source("R/view.R")
view_figure("outputs/part1/depth_profiles.png")
```

</details>

![Carbon accumulated down each core](../outputs/part1/cumulative_stock.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Re-read the workbook and rebuild Step 2 outputs.
source("scripts/01_field.R")
source("R/view.R")
view_figure("outputs/part1/cumulative_stock.png")
```

</details>

**Analyse:** Plot carbon concentration, density and cumulative stock using `ggplot2`.

**Read:** Each profile point represents a sampled slice; depth increases downward. Greater concentration alone does not guarantee greater stock: density and slice thickness matter too. The cumulative curve stops at sampled depth.

**Apply:** Identify changes down the profile and choose a consistent reporting depth.

[Profile and cumulative-curve methods](../R/field.R)

### Step 3 Compare other ground data

**Input:** Step 2's [core stocks](../outputs/part1/core_stocks.csv) and [layer stocks](../outputs/part1/layer_stocks.csv), plus [reference summaries](../data/reference/).

**Output:** Comparison figures and [community-to-reference table](../outputs/part1/comparison/community_vs_reference_groups.csv). **Used next:** These comparisons provide context; Part 2A returns to Step 2's core stocks.

![Comparison by collection](../outputs/part1/comparison/stock_by_dataset.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/stock_by_dataset.png")
```

</details>

![Comparison by ecosystem](../outputs/part1/comparison/stock_by_ecosystem.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/stock_by_ecosystem.png")
```

</details>

**Analyse:** Compare with precomputed summaries from NRCan peat, WoSIS, NPDB, Janousek and Riley. Match quantity, units, depth and depth reference.

**Read:** Thick segments show the middle half of reference values; dots show medians, diamonds show means, and dashed lines show your cores. Some ecosystem and study labels are unavailable. Riley contributes concentration only.

**Apply:** Put your results in context and choose relevant reference groups. These descriptive comparisons do not establish why groups differ or validate an illustrative core.

[Comparison methods](../R/compare.R) · [Provenance and source limits](../SOURCES.md)

<details>
<summary>More comparison figures</summary>

![Comparison by general region](../outputs/part1/comparison/stock_by_region.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/stock_by_region.png")
```

</details>

![Comparison by study or collection](../outputs/part1/comparison/stock_by_study.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/stock_by_study.png")
```

</details>

![Comparison by province or source region](../outputs/part1/comparison/stock_by_source_region.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/stock_by_source_region.png")
```

</details>

![Carbon concentration by collection](../outputs/part1/comparison/concentration_by_dataset.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/concentration_by_dataset.png")
```

</details>

![Bulk density by collection](../outputs/part1/comparison/density_by_dataset.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires core_stocks.csv and layer_stocks.csv from Step 2.
source("scripts/02_compare.R")
source("R/view.R")
view_figure("outputs/part1/comparison/density_by_dataset.png")
```

</details>

Study and source-region figures show the 20 largest groups by record count. Complete tables retain all groups, intervals and datums. Counts are profile records, not necessarily independent sites.

</details>

### Measurements and comparisons

Use the mean to summarise cores, sample SD to describe their spread, and `CV (%) = 100 × SD / mean` to express relative spread. R uses `mean()` and `stats::sd()`.

The two-core forest mean is 25.01 kg C/m²; the two-core wetland mean is 11.57. These describe this example, not the wider landscape. [View statistics](../outputs/part1/descriptive_statistics.csv).

### Data to retain

Keep core/sample IDs, ecosystem, slice boundaries, depth reference, density, carbon, units, measurement methods, coordinates/CRS and origin. Keep measured, estimated and illustrative values distinguishable. Retain QC flags alongside outputs.

## Part 2 Area based carbon stock analysis

### A Estimate totals by ecosystem

**Input:** Step 2's [core stocks and coordinates](../outputs/part1/core_stocks.csv), plus the [area boundary](../data/aoi/study_area.geojson).

**Output:** [Area totals](../outputs/part2/area/area_totals.csv), map figures, and [area checkpoint](../outputs/part2/area/area_analysis.rds). **Used next:** Part 2B reads that checkpoint, including the grid and four core stocks.

![Conceptual ecosystem areas](../outputs/part2/area/ecosystem_map.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires Step 2 core_stocks.csv and data/aoi/study_area.geojson.
source("scripts/03_area.R")
source("R/view.R")
view_figure("outputs/part2/area/ecosystem_map.png")
```

</details>

![Group mean stocks mapped across the area](../outputs/part2/area/group_stock_map.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires Step 2 core_stocks.csv and data/aoi/study_area.geojson.
source("scripts/03_area.R")
source("R/view.R")
view_figure("outputs/part2/area/group_stock_map.png")
```

</details>

**Analyse:** Apply each group's mean 0–30 cm stock to its area:

`Total (tonnes C) = mean stock (kg C/m²) × area (m²) / 1,000`

**Read:** `sf` calculates clipped cell areas; `terra` writes map grids. The example assigns ecosystem zones by the nearest core. They are hypothetical land cover.

**Apply:** Replace these zones with mapped ecosystems and representative sampling before estimating a real community inventory. [View totals](../outputs/part2/area/area_totals.csv).

[Area methods](../R/spatial.R)

### B Estimate between samples

**Input:** Part 2A's [area checkpoint](../outputs/part2/area/area_analysis.rds).

**Output:** The distance-weighted map and [spatial checkpoint](../outputs/part2/spatial/spatial_analysis.rds). This run also saves the separate simulation inspected in C. **Used next:** C inspects the simulation; Part 3 uses this distance-weighted map.

![Distance-weighted stock map](../outputs/part2/spatial/idw_stock_map.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires Part 2A area_analysis.rds; also creates the C simulation.
source("scripts/04_spatial.R")
source("R/view.R")
view_figure("outputs/part2/spatial/idw_stock_map.png")
```

</details>

**Analyse:** Give nearer cores more influence using inverse-distance weighting: `wᵢ = 1 / distanceᵢ²`; `stock(s) = Σ(wᵢ × stockᵢ) / Σwᵢ`. Base R calculates the weighted average; `sf::st_distance()` supplies distances.

**Read:** Smooth transitions show the weighting rule. Four cores do not establish spatial accuracy; this step does not fit kriging or a variogram.

**Apply:** Explore how sample placement affects the picture and where more measurements could help.

[Distance-weighting methods](../R/spatial.R)

### C Use landscape information

**Input:** The [spatial checkpoint](../outputs/part2/spatial/spatial_analysis.rds), figures and [holdout metrics](../outputs/part2/spatial/holdout_metrics.csv) already created in B.

**Output:** Inspect the saved model map and test results below. **Used next:** This is a separate simulation branch. Part 3 returns to B's distance-weighted map, rather than this regression map.

![Landscape model map](../outputs/part2/spatial/model_stock_map.png)

[Your figure caption]

<details>
<summary>Copy R code: inspect and view this figure</summary>

```r
# Inspect the checkpoint already created in B; do not rerun the model.
spatial_results <- readRDS("outputs/part2/spatial/spatial_analysis.rds")
stats::coef(spatial_results$final_model)
source("R/view.R")
view_figure("outputs/part2/spatial/model_stock_map.png")
```

</details>

![Held-out simulated and predicted stocks](../outputs/part2/spatial/observed_predicted.png)

[Your figure caption]

<details>
<summary>Copy R code: inspect and view this figure</summary>

```r
# Inspect the simulation test outputs already created in B.
read.csv("outputs/part2/spatial/holdout_metrics.csv")
source("R/view.R")
view_figure("outputs/part2/spatial/observed_predicted.png")
```

</details>

**Analyse:** Fit `stock ~ wetness + elevation` with `stats::lm()`. This separate exercise creates 160 simulated observations and simulated predictors: 112 train the model and 48 eastern observations test it.

**Read:** Points near the diagonal have smaller prediction errors. `RMSE = √mean((observed − predicted)²)` summarises error size. These are simulation results, not extra field observations or validation of real carbon stocks. The final demonstration map is refitted after testing.

**Apply:** With real data, select useful predictors and test predictions at locations excluded from fitting.

[Simulation, regression and spatial-holdout methods](../R/spatial.R)

## Part 3 Bayesian updating of prior maps

**Input:** Part 2B's distance-weighted values in the [spatial checkpoint](../outputs/part2/spatial/spatial_analysis.rds), plus an invented prior and assumed uncertainties.

**Output:** Updated maps, [summary](../outputs/part3/bayesian_summary.csv), and [cell values](../outputs/part3/bayesian_grid.csv). **Used next:** Inspect the change and uncertainty before replacing illustrative inputs.

![Prior, updated stock and change](../outputs/part3/bayesian_maps.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires Part 2B spatial_analysis.rds; updates its IDW map.
source("scripts/05_bayes.R")
source("R/view.R")
view_figure("outputs/part3/bayesian_maps.png")
```

</details>

![Assumed and updated uncertainty](../outputs/part3/uncertainty.png)

[Your figure caption]

<details>
<summary>Copy R code: make and view this figure</summary>

```r
# Requires Part 2B spatial_analysis.rds; updates its IDW map.
source("scripts/05_bayes.R")
source("R/view.R")
view_figure("outputs/part3/uncertainty.png")
```

</details>

**Analyse:** Combine a prior estimate and a new estimate, giving more weight to the estimate with less uncertainty. Base R implements the normal–normal update:

`Vpost = 1 / (1/Vprior + 1/Vnew)`

`Meanpost = Vpost × (Meanprior/Vprior + Meannew/Vnew)`

**Read:** `V` means variance, or SD². The example combines an invented prior map with the distance-weighted map. Prior SD = 5 and new-estimate SD = 3 kg C/m² are assumptions; updated SD is about 2.57.

**Apply:** Match quantities and use justified uncertainties before updating a real map. The formula assumes independent errors. Core-to-core SD is not automatically map uncertainty, and cell uncertainty does not give uncertainty for the whole-area total.

[Bayesian formula and assumptions](../R/spatial.R)

## 4 Case studies

### Fort Severn

[Your case-study text, figures and interpretation]

The real Fort Severn analysis has not been rerun in this illustrative course.

## Code and detailed methods

`readxl` imports the workbook; base R/`stats` calculate stocks, summaries, regression and updating; `sf` handles geometry and distance; `terra` writes rasters; `ggplot2` draws figures.

Run steps in order, or use `source("run_all.R")`. The [input/output guide](data-flow.md) shows each dependency. If you change the workbook, rerun Step 2 and the later steps; opening an existing image alone does not refresh it. [Field methods](../R/field.R), [comparison methods](../R/compare.R) and [spatial/Bayesian methods](../R/spatial.R) contain the full functions. Each spatial output folder includes `method_metadata.csv` with its formula, packages and assumptions.

The course rebuilds community comparisons and figures from [precomputed reference summaries](../SOURCES.md); it does not rebuild those summaries from raw reference profiles. Keep the original data, configuration and QC outputs with your results.
