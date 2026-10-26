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
