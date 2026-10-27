source("R/field.R")

expect_error <- function(expr) {
  failed <- inherits(tryCatch({force(expr);NULL},error=function(e)e),"error")
  if (!failed) stop("Expected an invalid input to be rejected.")
}
layers <- data.frame(core_id="unrelated-id",section_id=c("a","b"),
  top_cm=c(0,10),bottom_cm=c(10,30),bd_g_cm3=c(.1,.2),carbon_used_pct=c(50,10),
  ecosystem="Bog",group="Wetland cores")
metadata <- data.frame(core_id="unrelated-id",plot_id="p",ecosystem="Bog",group="Wetland cores",
  latitude=52,longitude=-84)
x <- calculate_core_stocks(layers,metadata)
# Independently: .1 g/cm3 * .5 * 10 cm * 10 = 5 kg/m2;
#              .2 g/cm3 * .1 * 20 cm * 10 = 4 kg/m2.
stopifnot(abs(x$cores$stock_0_30_kg_m2-9)<1e-10,
  abs(x$cores$stock_0_30_t_ha-90)<1e-10,
  is.na(x$cores$stock_0_100_kg_m2),x$cores$complete_sampled_profile)
zero <- layers; zero$carbon_used_pct[1] <- 0
stopifnot(calculate_core_stocks(zero,metadata)$cores$stock_0_30_kg_m2==4)
gap <- layers; gap$top_cm[2] <- 15
g <- calculate_core_stocks(gap,metadata)
stopifnot(is.na(g$cores$stock_0_30_kg_m2),!g$cores$complete_sampled_profile,
  abs(g$cores$sampled_stock_kg_m2-8)<1e-10)
deep <- layers; deep$bottom_cm[2] <- 150
d <- calculate_core_stocks(deep,metadata)
stopifnot(abs(d$cores$stock_0_100_kg_m2-23)<1e-10,
  abs(d$cores$stock_0_30_kg_m2-9)<1e-10)
overlap <- layers; overlap$top_cm[2] <- 9
expect_error(calculate_core_stocks(overlap,metadata))
missing <- layers; missing$bd_g_cm3[1] <- NA_real_
expect_error(calculate_core_stocks(missing,metadata))
invalid <- layers; invalid$carbon_used_pct[1] <- Inf
expect_error(calculate_core_stocks(invalid,metadata))
duplicate <- layers; duplicate$section_id[2] <- duplicate$section_id[1]
expect_error(calculate_core_stocks(duplicate,metadata))
expect_error(calculate_core_stocks(layers,rbind(metadata,metadata)))
# An example-specific correction must never change an unrelated user's row 7.
example_row <- data.frame(core_id="SIM-WET-01",section_id="SIM-WET-01-03",excel_row=7L,bd_g_cm3=.9)
stopifnot(validate_example_density_override(example_row,.09))
example_row$bd_g_cm3 <- .09
stopifnot(validate_example_density_override(example_row,.09))
own_row <- example_row; own_row$core_id <- "my-own-core"
expect_error(validate_example_density_override(own_row,.09))
stopifnot(validate_example_density_override(own_row,NA_real_))
wrong_section <- example_row; wrong_section$section_id <- "another-section"
expect_error(validate_example_density_override(wrong_section,.09))
wrong_density <- example_row; wrong_density$bd_g_cm3 <- .12
expect_error(validate_example_density_override(wrong_density,.09))
message("Field checks passed: units, zero carbon, fixed-depth coverage, gaps, overlap, missing inputs and guarded example correction.")

# A simple field table works with any core names and ecosystem labels.
flat <- data.frame(core_id=c("Ottawa-A","Ottawa-A","North-B","Prairie-C"),
  ecosystem=c("Marsh","Marsh","Tundra","Grassland"),latitude=c(45,45,68,51),
  longitude=c(-75,-75,-110,-106),top_cm=c(0,10,0,0),bottom_cm=c(10,30,30,20),
  bulk_density_g_cm3=c(.1,.2,.15,.3),organic_carbon_pct=c(50,10,20,5))
flat_path <- tempfile(fileext=".csv")
read_flat <- function(table, override=NA_real_) {
  utils::write.csv(table,flat_path,row.names=FALSE,na="")
  import_core_workbook(flat_path,density_override_J7=override)
}
flat_input <- read_flat(flat)
single <- read_flat(flat[3, , drop=FALSE])
stopifnot(nrow(single$layers)==1,nrow(single$cores)==1)
flat_stocks <- calculate_core_stocks(flat_input$layers,flat_input$cores)$cores
stopifnot(nrow(flat_input$layers)==4, nrow(flat_input$cores)==3,
  length(unique(flat_input$layers$section_id))==4,
  identical(sort(flat_input$cores$ecosystem),sort(c("Marsh","Tundra","Grassland"))),
  identical(sort(flat_input$cores$group),sort(c("Marsh","Tundra","Grassland"))),
  flat_stocks$stock_0_30_kg_m2[flat_stocks$core_id=="Ottawa-A"]==9,
  is.na(flat_stocks$stock_0_30_kg_m2[flat_stocks$core_id=="Prairie-C"]),
  identical(flat_input$depth_datum,"surface"))
expect_error(read_flat(flat[,names(flat)!="organic_carbon_pct"]))
expect_error(read_flat(flat,.09))
bad_metadata <- flat; bad_metadata$latitude[2] <- 46
expect_error(read_flat(bad_metadata))
bad_metadata <- flat; bad_metadata$ecosystem[2] <- "Swamp"
expect_error(read_flat(bad_metadata))
missing_coordinates <- flat
missing_coordinates$latitude[1:2] <- NA_real_
missing_coordinates$longitude[1:2] <- NA_real_
no_coords <- read_flat(missing_coordinates)
no_coord_stocks <- calculate_core_stocks(no_coords$layers,no_coords$cores)$cores
stopifnot(any(no_coords$qc$check=="missing_coordinates"),
  no_coord_stocks$stock_0_30_kg_m2[no_coord_stocks$core_id=="Ottawa-A"]==9)
bad_carbon <- flat; bad_carbon$organic_carbon_pct[1] <- NA_real_
expect_error(read_flat(bad_carbon))
bad_density <- flat; bad_density$bulk_density_g_cm3[1] <- "unknown"
expect_error(read_flat(bad_density))
optional_ids <- flat; optional_ids$sample_id <- paste0("sample-",1:4)
optional_ids$plot_id <- c("plot-A","plot-A","plot-B","plot-C")
with_ids <- read_flat(optional_ids)
stopifnot(identical(with_ids$layers$section_id,optional_ids$sample_id[order(optional_ids$core_id,optional_ids$top_cm)]))
optional_ids$sample_id[2] <- optional_ids$sample_id[1]
expect_error(read_flat(optional_ids))
unlink(flat_path)

# The simplified template reproduces all four corrected example stocks.
if (file.exists("data/input/field_data_template.csv")) {
  template <- import_core_workbook("data/input/field_data_template.csv")
  stopifnot(nrow(template$layers)==29,nrow(template$cores)==4)
  original <- import_core_workbook("data/input/community_cores.xlsm",.09)
  original_stocks <- calculate_core_stocks(original$layers,original$cores)$cores
  template_stocks <- calculate_core_stocks(template$layers,template$cores)$cores
  compare_columns <- setdiff(names(original_stocks),"plot_id")
  stopifnot(isTRUE(all.equal(original_stocks[compare_columns],template_stocks[compare_columns],tolerance=1e-10,
    check.attributes=FALSE)))
}
message("Simple field table checks passed: named columns, arbitrary IDs/habitats, core metadata, missing coordinates, required measurements and example reconciliation.")
