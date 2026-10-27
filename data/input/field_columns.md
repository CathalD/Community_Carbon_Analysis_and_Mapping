# Fill in your field data

Open `field_data_template.csv` in Excel or another spreadsheet program. It contains the 29 revised example slices and the confirmed density correction. Replace these rows with your own measurements. Keep the eight column names in row 1. Save as CSV, or as an Excel `.xlsx` file with this table on the first sheet.

Use one row for each measured soil slice. Repeat the core's name, ecosystem and coordinates on every row for that core. A core is one sampling location, with one or more slices down its depth.

| Column | What to enter |
| --- | --- |
| `core_id` | A unique name for each core, such as `Core-01`. Use the same name for all its slices. |
| `ecosystem` | Your ecosystem label, such as `Bog`, `Forest`, `Marsh`, `Tundra` or `Grassland`. Repeat exactly for each slice in that core. |
| `latitude` | WGS84 latitude in decimal degrees, such as `52.60`. Repeat for each slice. |
| `longitude` | WGS84 longitude in decimal degrees, such as `-84.40`. Locations west of Greenwich have a minus sign, including most of Canada. Repeat for each slice. |
| `top_cm` | Depth at the top of the slice, in centimetres below the ground surface. The first slice starts at `0`. |
| `bottom_cm` | Depth at the bottom of the slice, in centimetres. It must be greater than `top_cm`. |
| `bulk_density_g_cm3` | Dry bulk density in grams per cubic centimetre of the whole sampled volume. Enter `0.09`, for example. |
| `organic_carbon_pct` | Organic carbon as a percentage of dry sample mass. Enter 45% as `45`, not `0.45`. |

Measure depth from the ground surface, including the surface organic layer. Use this same starting point for every core. Each slice should begin where the previous one ends. The course calculates a 0–30 cm stock only when slices cover that full depth without gaps. It leaves missing deeper stocks blank.

These three rows show two cores. `Core-01` has two slices and `Core-02` has one:

```csv
core_id,ecosystem,latitude,longitude,top_cm,bottom_cm,bulk_density_g_cm3,organic_carbon_pct
Core-01,Marsh,52.60,-84.40,0,10,0.10,50
Core-01,Marsh,52.60,-84.40,10,30,0.20,10
Core-02,Grassland,52.61,-84.39,0,30,0.30,5
```

Use measured dry bulk density and organic carbon. This simple table does not convert loss on ignition into carbon. The stock calculation assumes the measurements describe the whole sample volume; correct your measurements first if coarse fragments need a separate adjustment.

Coordinates may be blank if you only want stocks and profile figures. A core needs both coordinates and a complete 0–30 cm stock to appear on the stock maps. Do not leave the example coordinates attached to your own measurements.

Optional columns `sample_id` and `plot_id` can retain your existing labels. Sample IDs must be unique across the file. Otherwise the code creates them. A plot ID must stay the same on every row for its core.

Save your table as `data/input/my_field_data.csv`, then run `source("run_my_data.R")` from the course R project. This route automatically disables the older workbook's J7 correction. If you save an Excel `.xlsx` file instead, put the table on its first sheet and change only the `my_field_data` path at the top of `run_my_data.R`.
