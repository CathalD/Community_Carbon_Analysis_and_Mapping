# Start here

**You need:** [R](https://cran.r-project.org/) and [RStudio Desktop](https://docs.posit.co/ide/user/download.html). Install R first, then RStudio. Run the example before adding your own soil measurements and study boundary.

The route is:

**Field sheet → checked samples → carbon stocks and profiles → sample locations → a simple carbon map.**

## 1 Download and open the course

Paste this into the RStudio **Console**:

```r
install.packages("usethis") # Once only.
usethis::use_course(
  "https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip"
)
```

Open `Community_Carbon_Analysis_and_Mapping.Rproj` in the downloaded folder. Keep the folder structure unchanged.

## 2 Run the supplied example

```r
source("run_example.R")
```

`source()` runs a saved R script. Setup automatically installs missing packages: `readxl`, `ggplot2`, `sf` and `terra`. The first run needs internet access. Wait for installation to finish.

The example uses four illustrative cores in western James Bay. Open `outputs/example/` in your file browser, or paste:

```r
source("R/view.R")
view_figure("outputs/example/part1/depth_profiles.png")
```

The [workshop](workshop.md) shows all five example figures with their code.

## 3 Put your measurements in the field sheet

Copy [field_data_template.csv](../data/input/field_data_template.csv). Save the copy as **`my_field_data.csv`** inside **`data/input/`**.

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

Start depth at **0 at the ground surface**, including the organic surface layer. Slices must not overlap. A 0–30 cm stock needs complete coverage from 0 to 30 cm; a gap or a shallower core gives a missing value for that depth.

Keep numbers numeric: enter `10`, not `10 cm`. Save as a normal comma-separated **CSV** file. This simple route does not need macros or extra worksheets.

## 4 Put your study boundary in the AOI folder

**AOI** means **area of interest**: the boundary of the place you want to map.

Save your polygon boundary as **`my_aoi.geojson`** inside **`data/aoi/`**. GeoJSON stores the outline and its geographic position. Export an existing boundary as GeoJSON with a known coordinate system.

To draw a boundary, open [geojson.io](https://geojson.io/) in your browser. Zoom to your study area, choose the polygon tool, click the corners, and close the outline. Download/export GeoJSON, rename it `my_aoi.geojson`, and move it to `data/aoi/`. See the [tool documentation](https://github.com/mapbox/geojson.io).

Use a local area containing your samples. Changing the boundary does **not** move the example cores: replace their coordinates too. The map needs at least **two different locations inside the boundary**, each with a complete 0–30 cm stock. A large area with few samples is poorly supported.

## 5 Run your data and inspect the results

```r
source("run_my_data.R")
```

With these filenames, no code changes are needed. Each run refreshes `outputs/my_data/` and preserves your field sheet.

| Open this file | Check or use it for |
| --- | --- |
| `part1/qc_flags.csv` | Data checks to review first. |
| `part1/core_stocks.csv` | Sampled stock and equal-depth stocks for each core. |
| `part1/descriptive_statistics.csv` | Group means and spread. |
| `part1/common_depth_stock.png` | Compare cores over 0–30 cm. |
| `part1/depth_profiles.png` | Read density and carbon down each core. |
| `part1/cumulative_stock.png` | See carbon added with depth. |
| `maps/sample_map.png` | Check where the usable samples lie. |
| `maps/idw_stock_map.png` | Explore estimates between the samples. |
| `maps/idw_stock_0_30.tif` | Open the map in mapping software. |
| `maps/excluded_cores.csv` | See samples left out of mapping and why. |
| `maps/mapping_metadata.csv` | Check mapping settings. |

These paths are inside `outputs/my_data/`. Map shading is an estimate, not an extra measurement.

After editing inputs, rerun `source("run_my_data.R")` to refresh figures.

## If something stops

- **File not found:** open the course `.Rproj`; check the exact filenames and folders above.
- **Missing column or invalid number:** compare the field sheet with the eight headers and units above.
- **No mapped samples:** check negative longitude, your boundary, and complete 0–30 cm coverage.
- **Package installation failed:** read the installation message. Some systems need extra system libraries for `sf` or `terra`; ask your instructor or computer support for help.

Finish these five figures before trying the [optional full course](workshop.md#optional-full-illustrative-course).
