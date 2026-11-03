# Harmonise profile layers to the map depth intervals for comparison.
#
# Method (agreed 15 Oct 2026): thickness-weighted mean of the layers that
# overlap each interval, computed only where the profile's layers cover the
# whole interval with no gaps and no overlaps. No extrapolation, no splines.
#   * Concentrations (SOC g/kg, BD g/cm3): overlap-weighted mean.
#   * Stocks (kg C/m2) for 0-30 and 0-100 cm: sum of each layer's own stock
#     value times the fraction of the layer inside the interval, using the
#     dataset's own stock fields only (no new bulk-density gap filling).
#   * Two depth datums, kept separate (they are not interchangeable):
#       surface - depth from the top of the profile, organic surface layers
#                 included (NPDB shifted so its profile top is 0; others as
#                 prepared). Used with organic-inclusive stock maps.
#       mineral - organic surface layers (WoSIS organic_surface = TRUE; NPDB
#                 horizons above 0 cm, i.e. LFH) removed and the first remaining
#                 layer set to 0. WoSIS and NPDB only. Used with mineral-soil
#                 products (SoilGrids removes overlying organic layers).
#   * Riley bulk density units are unverified: excluded from quantitative use
#     (kept as bd_unverified in the harmonised layers).
#   * Janousek profiles whose recorded thickness disagrees with their
#     boundaries (> 0.01 cm) are excluded from depth integration.
#
# Inputs : prepared layer CSVs in Carbon_Resources/ (parent of the repo)
# Output : ../../data/derived/profile_intervals.csv  (long: profile x interval x variable x datum)
# Run    : /Library/Frameworks/R.framework/Resources/bin/Rscript 01_profile_intervals.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
local_dir <- normalizePath(file.path(here, "..", "..", "..", ".."))
out_dir <- normalizePath(file.path(here, "..", "..", "data", "derived"), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) read_csv(file.path(local_dir, f), show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

# ---------------------------------------------------------------- common units
# soc_gkg: organic carbon g/kg; bd: g/cm3; stock_kgm2: layer carbon stock kg C/m2
wosis <- rd("wosis_canada_layers_v2.csv") |>
  transmute(dataset = "WoSIS", profile_uid, top = upper_depth, bot = lower_depth,
            organic_surface = as.logical(organic_surface),
            soc_gkg = OrgC, bd = BDOD, stock_kgm2 = OrgC_stock_kgm2)
npdb <- rd("npdb_canada_layers_v1_upload_fixed.csv") |>
  transmute(dataset = "NPDB", profile_uid, top = top_depth_cm, bot = bot_depth_cm,
            organic_surface = top_depth_cm < 0,   # LFH horizons sit above the mineral surface (0 cm)
            soc_gkg = carbon_pct * 10, bd = bulk_density, stock_kgm2 = stock_Mgha / 10)
jan_raw <- rd("janousek_layers_v1.csv")
jan_bad <- jan_raw |> filter(abs(layer_thickness_cm - (lower_depth - upper_depth)) > 0.01) |> distinct(profile_uid)
message(sprintf("Janousek profiles excluded for thickness/boundary mismatch: %d", nrow(jan_bad)))
others <- bind_rows(
  rd("canpeat_layers_v1.csv") |>
    transmute(dataset = "NRCan peat", profile_uid, top = upper_depth, bot = lower_depth,
              soc_gkg = OrgC * 10, bd = BDOD,
              # provider sample carbon mass (SAMP_CARB_MGHA, Mg C/ha) -> kg C/m2
              stock_kgm2 = TOTC_Stock / 10),
  jan_raw |> anti_join(jan_bad, by = "profile_uid") |>
    transmute(dataset = "Janousek", profile_uid, top = upper_depth, bot = lower_depth,
              soc_gkg = OrgC_pct * 10, bd = BDOD, stock_kgm2 = OrgC_Stock_kgm2),
  rd("riley_layers_v1.csv") |>
    transmute(dataset = "Riley", profile_uid, top = upper_depth_cm, bot = lower_depth_cm,
              soc_gkg = organic_carbon_pct * 10,
              bd = NA_real_, bd_unverified = dry_bulk_density_source,   # units unverified: not used
              stock_kgm2 = NA_real_)
)

surface <- bind_rows(
  wosis |> select(-organic_surface),
  npdb |> group_by(profile_uid) |> mutate(shift = min(top, na.rm = TRUE)) |> ungroup() |>
    mutate(top = top - shift, bot = bot - shift) |> select(-shift, -organic_surface),
  others
) |> mutate(datum = "surface")

mineral <- bind_rows(wosis, npdb) |>
  filter(!organic_surface) |>
  group_by(profile_uid) |> mutate(shift = min(top, na.rm = TRUE)) |> ungroup() |>
  mutate(top = top - shift, bot = bot - shift) |> select(-shift, -organic_surface) |>
  mutate(datum = "mineral")

layers <- bind_rows(surface, mineral) |>
  filter(!is.na(top), !is.na(bot), bot > top)

# The harmonised layer table is an in-memory intermediate; a 2026-10-16 snapshot is kept in
# archive/derived-snapshots/ for reference.

# ---------------------------------------------------------------- interval logic
intervals <- tribble(
  ~interval, ~itop, ~ibot,
  "0-5", 0, 5, "5-15", 5, 15, "15-30", 15, 30, "30-60", 30, 60, "60-100", 60, 100,
  "0-30", 0, 30, "0-100", 0, 100
)

# For one profile's layers and one interval: coverage check, then aggregate.
summarise_interval <- function(L, itop, ibot) {
  ov <- pmax(0, pmin(L$bot, ibot) - pmax(L$top, itop))
  keep <- ov > 0
  if (!any(keep)) return(NULL)
  L <- L[keep, ]; ov <- ov[keep]
  # Full, non-overlapping coverage: overlaps sum to the interval thickness and
  # the union of clipped layers spans it.
  span <- ibot - itop
  clipped_top <- pmax(L$top, itop); clipped_bot <- pmin(L$bot, ibot)
  o <- order(clipped_top)
  gaps_or_overlaps <- any(abs(clipped_top[o][-1] - clipped_bot[o][-length(o)]) > 0.5)
  covered <- abs(sum(ov) - span) <= 0.5 && min(clipped_top) <= itop + 0.5 &&
    max(clipped_bot) >= ibot - 0.5 && !gaps_or_overlaps
  if (!covered) return(NULL)
  frac <- ov / (L$bot - L$top)
  tibble(
    soc_gkg = if (all(!is.na(L$soc_gkg))) sum(L$soc_gkg * ov) / span else NA_real_,
    bd = if (all(!is.na(L$bd))) sum(L$bd * ov) / span else NA_real_,
    stock_kgm2 = if (all(!is.na(L$stock_kgm2))) sum(L$stock_kgm2 * frac) else NA_real_,
    n_layers = nrow(L)
  )
}

res <- layers |>
  group_by(datum, dataset, profile_uid) |>
  group_modify(function(L, key) {
    map_dfr(seq_len(nrow(intervals)), function(i) {
      r <- summarise_interval(L, intervals$itop[i], intervals$ibot[i])
      if (is.null(r)) NULL else mutate(r, interval = intervals$interval[i])
    })
  }) |>
  ungroup()

long <- res |>
  pivot_longer(c(soc_gkg, bd, stock_kgm2), names_to = "variable", values_to = "observed") |>
  filter(!is.na(observed)) |>
  # concentrations are compared at the five map intervals; stocks at 0-30 / 0-100
  filter((variable %in% c("soc_gkg", "bd") & interval %in% c("0-5", "5-15", "15-30", "30-60", "60-100")) |
           (variable == "stock_kgm2" & interval %in% c("0-30", "0-100")))

write_csv(long, file.path(out_dir, "profile_intervals.csv"))
message(sprintf("layers: %s | interval values: %s", format(nrow(layers), big.mark = ","), format(nrow(long), big.mark = ",")))
print(long |> count(datum, dataset, variable, interval) |> pivot_wider(names_from = interval, values_from = n), n = 50, width = 200)
