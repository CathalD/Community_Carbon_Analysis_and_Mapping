# Run the full course first, then these focused mathematical and data checks.
source("tests/check_field.R")
source("tests/check_spatial.R")
core <- read.csv("outputs/part1/core_stocks.csv")
own <- read.csv("outputs/part1/comparison/community_matched_depth_values.csv")
stock <- subset(own, variable == "stock_kgm2" & interval == "0-30" & eligible)
stopifnot(nrow(stock) == nrow(core),
  max(abs(stock$value - core$stock_0_30_kg_m2[match(stock$core_id,core$core_id)])) < 1e-9)
depth100 <- subset(own, variable == "stock_kgm2" & interval == "0-100")
stopifnot(all(!depth100$eligible[depth100$core_id %in% core$core_id[!core$complete_to_100cm]]))
expected <- c("depth_profiles.png", "common_depth_stock.png", "cumulative_stock.png")
stopifnot(all(file.exists(file.path("outputs/part1", expected))))
message("Course checks passed.")
