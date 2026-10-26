source("config.R")
source("R/field.R")

core_input <- import_core_workbook(course_config$workbook_path,
  density_override_J7=course_config$density_override_J7)
field <- analyse_core_data(core_input,out="outputs/part1",data_origin=course_config$data_origin)
print(field$cores[c("core_id","ecosystem","stock_0_30_kg_m2",
  "stock_0_100_kg_m2","sampled_stock_kg_m2")],row.names=FALSE)
message("Field outputs saved to outputs/part1. Review qc_flags.csv alongside the figures.")
