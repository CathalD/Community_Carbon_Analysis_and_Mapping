# Your area boundary

An **area of interest (AOI)** is the boundary of the area you want to map. It must be a closed shape, such as a community study area or watershed. Points and lines are not boundaries.

1. Obtain or draw your boundary in a GIS, then export it as **GeoJSON**. GeoJSON is a map file containing the boundary coordinates. Choose **WGS84 / EPSG:4326** when exporting; this stores longitude and latitude in degrees. Other correctly identified coordinate systems are also accepted.
2. Save the file in this folder as **`my_aoi.geojson`**. Keep the `.geojson` extension. A boundary with several separate polygons is accepted; they are treated as one study area.
3. The supplied **`run_my_data.R`** already looks for `data/aoi/my_aoi.geojson`. If you choose another file name, edit the `my_area` line at the top of that script. The `my_field_data` line directly above it names your field data file.
4. Enter your own core locations in the [field table](../input/field_columns.md) in **WGS84 decimal degrees**: latitude is positive and longitude is negative in Canada. The mapping step needs at least two cores at different locations inside the boundary, each with complete measurements from 0 to 30 cm below the ground surface. Save the completed table as **`data/input/my_field_data.csv`**.
5. Run `source("run_my_data.R")`. Open **`outputs/my_data/maps/sample_map.png`** first to check that the boundary and core locations are correct. Open **`idw_stock_map.png`** in the same folder to see the simple carbon estimate between locations.

Changing the boundary does not relocate the example cores. To map another area of Canada, use a boundary and core coordinates from that same area. Missing coordinates, incomplete 0–30 cm stocks and points outside the boundary are listed in **`excluded_cores.csv`**; their field calculations remain available.

The grid size is chosen automatically so small and large areas can run on an ordinary computer. Larger areas have coarser map cells. The chosen size is recorded in **`mapping_metadata.csv`**. The map estimates values from distance to the cores; it does not establish the map's accuracy or produce a verified carbon inventory.

`study_area.geojson` is the supplied Western James Bay example. Keep it for the supplied example course; use `my_aoi.geojson` for your own boundary.
