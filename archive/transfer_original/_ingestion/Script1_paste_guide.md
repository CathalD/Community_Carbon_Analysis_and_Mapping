# Script 1 paste instructions

The source scripts have not been edited. Run Script 1 from the Transfer_Oct29 folder.

## 1. Replace the CuOrgC_plot block

Replace from `CuOrgC_plot = ggplot(...)` through the line that displays `CuOrgC_plot`.

```r
CuOrgC_plot <- ggplot(sample_sheet, aes(
  x = `Bottom depth (cm)`,
  y = `Cumulative carbon to base (kg C/m²)`,
  colour = `Core ID`
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(title = "Cumulative carbon stock by depth",
       x = "Depth (cm)", y = "Cumulative carbon stock (kg C/m²)") +
  theme_minimal()

CuOrgC_plot
```

## 2. Replace the core_summary block

Replace from `core_summary <- sample_sheet` through `print(core_summary)` with the following. A 30–60 cm section contributes its stock / 30 × 20 to the 0–50 group, added to the stock above 30 cm. The group cutoffs are: 0–30 requires reaching 30; 0–50 requires deeper than 30; 0–100 requires deeper than 50. Total retains sampled profile stock.

```r
# Replace the existing core_summary block with this block.

stock_to_depth <- function(top, bottom, stock, limit) {
  if (any(!is.finite(top) | !is.finite(bottom) | bottom <= top)) return(NA_real_)
  maximum <- max(bottom)

  if (limit == 30 && maximum < 30) return(NA_real_)
  if (limit == 50 && maximum <= 30) return(NA_real_)
  if (limit == 100 && maximum <= 50) return(NA_real_)

  end <- min(limit, maximum)
  centimetres <- pmax(0, pmin(bottom, end) - pmax(top, 0))
  keep <- centimetres > 0
  if (!any(keep) || any(!is.finite(stock[keep]) | stock[keep] < 0)) return(NA_real_)

  # Check that the measured sections join without gaps or overlaps.
  a <- pmax(top[keep], 0)
  b <- pmin(bottom[keep], end)
  o <- order(a)
  if (abs(min(a)) > 1e-8 || abs(max(b) - end) > 1e-8 ||
      any(abs(a[o][-1] - head(b[o], -1)) > 1e-8)) return(NA_real_)

  carbon_per_cm <- stock[keep] / (bottom[keep] - top[keep])
  sum(carbon_per_cm * centimetres[keep])
}

core_summary <- sample_sheet %>%
  group_by(`Core ID`) %>%
  summarise(
    `Total Sections` = n(),
    `Max Depth (cm)` = max(`Bottom depth (cm)`),
    `Avg Bulk Density` = mean(`Bulk density (g/cm³)`, na.rm = TRUE),
    `Avg Carbon (%)` = mean(`Organic carbon (%)`, na.rm = TRUE),
    `C Stock 0-30cm` = stock_to_depth(`Top depth (cm)`, `Bottom depth (cm)`,
                                     `Stock, full thickness (kg C/m²)`, 30),
    `C Stock 0-50cm` = stock_to_depth(`Top depth (cm)`, `Bottom depth (cm)`,
                                     `Stock, full thickness (kg C/m²)`, 50),
    `C Stock 0-100cm` = stock_to_depth(`Top depth (cm)`, `Bottom depth (cm)`,
                                      `Stock, full thickness (kg C/m²)`, 100),
    `Total Carbon Stock` = stock_to_depth(`Top depth (cm)`, `Bottom depth (cm)`,
                                         `Stock, full thickness (kg C/m²)`, Inf),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(core_summary)
```

## 3. Fill the empty Part 3 section

Paste after core_summary is created and before the mapping export section. The data/reference folder is already prepared locally.

```r
# Paste into the empty Part 3 section, after core_summary has been created.
library(readr)

reference <- read_csv("data/reference/reference_depth_values.csv.gz",
                      show_col_types = FALSE) %>%
  filter(variable == "stock_kgm2", datum == "surface")

field <- bind_rows(
  core_summary %>% transmute(core_id = `Core ID`, interval = "0-30", value = `C Stock 0-30cm`),
  core_summary %>% transmute(core_id = `Core ID`, interval = "0-50", value = `C Stock 0-50cm`),
  core_summary %>% transmute(core_id = `Core ID`, interval = "0-100", value = `C Stock 0-100cm`),
  core_summary %>% transmute(core_id = `Core ID`, interval = "total", value = `Total Carbon Stock`)
) %>%
  filter(is.finite(value))

depth_order <- c("0-30", "0-50", "0-100", "total")
reference$interval <- factor(reference$interval, levels = depth_order)
field$interval <- factor(field$interval, levels = depth_order)

context_plot <- ggplot(reference, aes(x = dataset, y = value)) +
  geom_boxplot(fill = "grey90", outlier.alpha = 0.15) +
  geom_hline(data = field, aes(yintercept = value, colour = core_id),
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ interval, scales = "free_y") +
  labs(title = "Community cores and reference carbon stocks",
       x = NULL, y = "Carbon stock (kg C/m²)", colour = "Community core",
       caption = "Stocks are capped at each depth. Shorter eligible cores retain their measured stock.\nTotal represents the sampled profile; sampled depths vary.") +
  theme_minimal()

print(context_plot)

dir.create("outputs/script1", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/script1/carbon_stock_context.png", context_plot,
       width = 11, height = 8, dpi = 200)
write_csv(field, "outputs/script1/community_depth_groups.csv")
```

The reference bundle has the five previously used collections, source core/layer files, the four agreed depth groups, summary tables and resources.csv. Riley contributes carbon concentration only because its density units are unverified. The bundle is local staging; no GitHub upload has been performed. Source terms are retained in resources.csv and README.md.

For GitHub, the intended folder is Transfer_OCT_29/data/reference on branch Transfer_Oct. Keep this folder structure when downloading, then run Script 1 with Transfer_OCT_29 as the working directory.
