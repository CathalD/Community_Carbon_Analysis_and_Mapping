# Community Carbon Stock Analysis and Mapping

A university workshop that follows soil carbon data from field measurements to spatial estimates and Bayesian map updating. Work through the complete analysis using a supplied illustrative example from western James Bay, then practise the field calculations and mapping with your own data.

[Start here: files and commands](docs/start_here.md) · [Workshop: explanations, figures and code](docs/workshop.md) · [Explainer Google Doc](https://docs.google.com/document/d/1qflFw6CH69feNRr87Oi1SbHl0sY5bcXb-TuFI6KUtaU/edit)

## Download and run the workshop

In R or RStudio, paste:

```r
install.packages("usethis") # Once only.
usethis::use_course(
  "https://github.com/CathalD/Community_Carbon_Analysis_and_Mapping/archive/refs/heads/workshop-mvp.zip"
)
```

Open `Community_Carbon_Analysis_and_Mapping.Rproj` in the downloaded folder. In the **Console**, run:

```r
source("run_all.R")
```

`source()` runs a saved R script. This command checks packages and runs every workshop analysis. The first run needs internet access to install missing packages. The [start guide](docs/start_here.md) also gives the commands for working through one stage at a time.

## Workshop sequence

| Stage | Analysis | Results |
| --- | --- | --- |
| Part 1 | Check field data; calculate stocks, depth profiles and descriptive statistics; compare with reference collections, studies, ecosystems and regions. | `outputs/part1/` |
| Part 2A | Apply ecosystem means to areas and calculate totals. | `outputs/part2/area/` |
| Parts 2B–C | Map stocks with inverse-distance weighting; fit a landscape model and assess predictions using held-out data. | `outputs/part2/spatial/` |
| Part 3 | Combine a prior map and a new estimate using Bayesian updating. | `outputs/part3/` |

Each section links its inputs, runnable code, example figures and interpretation. Ecosystem areas, landscape data and the prior are illustrative teaching inputs.

## Practise with your data

Save measurements as `data/input/my_field_data.csv` and a matching boundary as `data/aoi/my_aoi.geojson`, then run `source("run_my_data.R")`. This separate practice script covers field stocks, profiles and IDW; it saves to `outputs/my_data/`. The [guide](docs/start_here.md) explains the columns, units and additional inputs needed to adapt later analyses.

Code is [MIT licensed](LICENSE); see [source details](SOURCES.md) for third-party data.
