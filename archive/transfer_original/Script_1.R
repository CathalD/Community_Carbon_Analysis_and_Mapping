
#===============================================================

# -------- PArt 1 - Ingest the dataset and organize the results

#===============================================================

#If packages not installed - please do so once:

#install.packages("janitor")
#install.packages("ggplot2")
#install.packages("readxl")

# ---- Step 1 ------

# ===================


#Configuring the data:

# If importing as the complete excel file run this:

#Change to the pathname of the workbook:
library(readxl)
workbook_path = "data/community_carbon_outputs_updated.xlsm"

read_excel(workbook_path)
excel_sheets(workbook_path)
#Import data for each core collected
core_sheet = read_excel(workbook_path, sheet = "2. Core Log", skip = 2)
head(core_sheet)
tail(core_sheet)
#Import the indivdual samples
sample_sheet <- read_excel(workbook_path, sheet = "3. Peat Data", skip = 2)
head(sample_sheet)
tail(sample_sheet)

#Clean up datasheet -
library(janitor)

# This strips out rows that contain nothing but NA values
core_sheet   <- core_sheet %>% remove_empty("rows")

library(dplyr)
# Keep rows only if they have a valid Plot ID and Core ID
sample_sheet <- sample_sheet %>%
  filter(!is.na(`Plot ID`) & !is.na(`Core ID`))



#===============================================================

# ------- Step 2 --------

#========================

# summarize the data

# Lets start with non-spatial data
# Lets look at the each core and compare depths, bulk density, carbon%
# First define variables:

OrgC = sample_sheet$`Organic carbon (%)`
BD = sample_sheet$`Bulk density (g/cm³)`
sample_depth = ((sample_sheet$`Bottom depth (cm)`)+(sample_sheet$`Top depth (cm)`))/2
CuOrgC = sample_sheet$`Cumulative carbon to base (kg C/m²)`
#First visually:
library(ggplot2)
OrgC_Plot <- ggplot(sample_sheet, aes(x = OrgC, y = sample_depth, group = `Core ID`, color = `Core ID`)) +
  # Add lines connecting the depths for each core
  geom_path(linewidth = 1) + 
  # Add points for each individual sample measurement
  geom_point(size = 2.5) +
  # REVERSE THE Y-AXIS (0 at the top, deepest at the bottom)
  scale_y_reverse(breaks = seq(0, 100, by = 10)) + 
  # Put the x-axis at the top if you prefer standard soil profile looks (optional)
  scale_x_continuous(position = "top") + 
  # Clean up labels and apply a minimalist theme
  labs(
    title = "Peat Organic Carbon Profile by Depth",
    x = "Organic Carbon (%)",
    y = "Depth (cm)"
  ) +
  theme_minimal()
BD_plot <- ggplot(sample_sheet, aes(x = BD, y = sample_depth, group = `Core ID`, color = `Core ID`)) +
  # Add lines connecting the depths for each core
  geom_path(linewidth = 1) + 
  # Add points for each individual sample measurement
  geom_point(size = 2.5) +
  # REVERSE THE Y-AXIS (0 at the top, deepest at the bottom)
  scale_y_reverse(breaks = seq(0, 100, by = 10)) + 
  # Put the x-axis at the top if you prefer standard soil profile looks (optional)
  scale_x_continuous(position = "top") + 
  # Clean up labels and apply a minimalist theme
  labs(
    title = "Peat Organic Carbon Profile by Depth",
    x = "Bulk Density (g/cm^2)",
    y = "Depth (cm)"
  ) +
  theme_minimal()
OrgC_Plot
BD_plot

# Lastly, show carbon accumulation curves for each:
CuOrgC_plot = ggplot(sample_sheet, aes(x = CuOrgC, y = sample_depth, group = `Core ID`, color = `Core ID`)) +
  # Add lines connecting the depths for each core
  geom_path(linewidth = 1) + 
  # Add points for each individual sample measurement
  geom_point(size = 2.5) +
  # REVERSE THE Y-AXIS (0 at the top, deepest at the bottom)
  scale_y_reverse(breaks = seq(0, 100, by = 10)) + 
  # Put the x-axis at the top if you prefer standard soil profile looks (optional)
  scale_x_continuous(position = "top") + 
  # Clean up labels and apply a minimalist theme
  labs(
    title = "Peat Organic Carbon Profile by Depth",
    x = "Cumulative Organic Carbon (kg/m^2)",
    y = "Depth (cm)"
  ) +
  theme_minimal()

#Optionl vertical version

CuOrgC_plot <- ggplot(sample_sheet, aes(
  x = `Bottom depth (cm)`,
  y = `Cumulative carbon to base (kg C/m²)`,
  colour = `Core ID`
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "Cumulative carbon stock by depth",
    x = "Depth (cm)",
    y = "Cumulative carbon stock (kg C/m²)"
  ) +
  theme_minimal()

CuOrgC_plot

# And some summary stats to go with it:
# Calculate summary stats per Core ID
library(dplyr)

# Calculate summary stats per Core ID with depth-specific carbon stocks

stock_to_depth <- function(top, bottom, stock, limit) {
  if (any(!is.finite(top) | !is.finite(bottom) | bottom <= top)) {
    return(NA_real_)
  }
  
  maximum <- max(bottom)
  
  if (limit == 30 && maximum < 30) return(NA_real_)
  if (limit == 50 && maximum <= 30) return(NA_real_)
  if (limit == 100 && maximum <= 50) return(NA_real_)
  
  end <- min(limit, maximum)
  centimetres <- pmax(0, pmin(bottom, end) - pmax(top, 0))
  keep <- centimetres > 0
  
  if (!any(keep) ||
      any(!is.finite(stock[keep]) | stock[keep] < 0)) {
    return(NA_real_)
  }
  
  # Check for gaps or overlapping sections.
  a <- pmax(top[keep], 0)
  b <- pmin(bottom[keep], end)
  o <- order(a)
  
  if (abs(min(a)) > 1e-8 ||
      abs(max(b) - end) > 1e-8 ||
      any(abs(a[o][-1] - head(b[o], -1)) > 1e-8)) {
    return(NA_real_)
  }
  
  carbon_per_cm <- stock[keep] / (bottom[keep] - top[keep])
  sum(carbon_per_cm * centimetres[keep])
}



core_summary <- sample_sheet %>%
  group_by(`Core ID`) %>%
  summarise(
    `Total Sections`      = n(),
    `Max Depth (cm)`      = max(`Bottom depth (cm)`, na.rm = TRUE),
    `Avg Bulk Density`    = mean(`Bulk density (g/cm³)`, na.rm = TRUE),
    `Avg Carbon (%)`      = mean(`Organic carbon (%)`, na.rm = TRUE),
    
    # Carbon Stocks at specific depth horizons
    `C Stock 0-30cm` = stock_to_depth(
      `Top depth (cm)`, `Bottom depth (cm)`,
      `Stock, full thickness (kg C/m²)`, 30
    ),
    
    `C Stock 0-50cm` = stock_to_depth(
      `Top depth (cm)`, `Bottom depth (cm)`,
      `Stock, full thickness (kg C/m²)`, 50
    ),
    
    `C Stock 0-100cm` = stock_to_depth(
      `Top depth (cm)`, `Bottom depth (cm)`,
      `Stock, full thickness (kg C/m²)`, 100
    ),
    
    `Total Carbon Stock` = stock_to_depth(
      `Top depth (cm)`, `Bottom depth (cm)`,
      `Stock, full thickness (kg C/m²)`, Inf
    )
  ) %>%
  # Round all numeric metrics to 2 decimal places
  mutate(across(where(is.numeric), ~ round(., 2)))

# View the new summary table with depth metrics
print(core_summary)


#=============================

# --- Export Results for Part 2

# ============================
library(dplyr)
library(readr)

# 1. Merge the coordinates and IDs from core_sheet into your summary
mapping_data <- core_summary %>%
  left_join(core_sheet, by = "Core ID") %>%
  # Relocate the spatial and identifier columns to the very front for clean mapping setup
  relocate(`Core ID`, `Plot ID`, Latitude, Longitude)

# 2. Export the final data frame to a mapping-ready CSV file
write_csv(mapping_data, "peat_cores_for_mapping.csv")

# Print a preview to ensure everything merged cleanly
head(mapping_data)

if (!dir.exists("Draft1/data")) {
  dir.create("Draft1/data", recursive = TRUE)
}

# 2. Export the final mapping dataset directly into your data folder
write_csv(mapping_data, "Draft1/data/mapping_data.csv")




# ========================================================

# ---------- Step 3 ----------------

# ==================================

# Bring in the external data

library(readr)

reference <- read_csv(
  "data/reference/reference_depth_values.csv.gz",
  show_col_types = FALSE
) %>%
  filter(variable == "stock_kgm2", datum == "surface")

field <- bind_rows(
  core_summary %>% transmute(
    core_id = `Core ID`, interval = "0-30",
    value = `C Stock 0-30cm`
  ),
  core_summary %>% transmute(
    core_id = `Core ID`, interval = "0-50",
    value = `C Stock 0-50cm`
  ),
  core_summary %>% transmute(
    core_id = `Core ID`, interval = "0-100",
    value = `C Stock 0-100cm`
  ),
  core_summary %>% transmute(
    core_id = `Core ID`, interval = "total",
    value = `Total Carbon Stock`
  )
) %>%
  filter(is.finite(value))

depth_order <- c("0-30", "0-50", "0-100", "total")
reference$interval <- factor(reference$interval, levels = depth_order)
field$interval <- factor(field$interval, levels = depth_order)

context_plot <- ggplot(reference, aes(x = dataset, y = value)) +
  geom_boxplot(fill = "grey90", outlier.alpha = 0.15) +
  geom_hline(
    data = field,
    aes(yintercept = value, colour = core_id),
    linetype = "dashed", linewidth = 0.6
  ) +
  facet_wrap(~ interval, scales = "free_y") +
  labs(
    title = "Field cores and reference carbon stocks",
    x = NULL, y = "Carbon stock (kg C/m²)",
    colour = "Community core"
  ) +
  theme_minimal()

print(context_plot)

dir.create("outputs/script1", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "outputs/script1/carbon_stock_context.png",
  context_plot, width = 11, height = 8, dpi = 200
)

write_csv(field, "outputs/script1/community_depth_groups.csv")

ecosystem_stock <- readr::read_csv(
  "data/reference/ecosystem_carbon_stocks.csv",
  show_col_types = FALSE
) %>%
  filter(ecosystem != "Not recorded") %>%
  mutate(interval = factor(
    interval,
    levels = c("0-30", "0-50", "0-100", "total")
  ))

ecosystem_counts <- ecosystem_stock %>%
  count(ecosystem, interval)

ecosystem_plot <- ggplot(
  ecosystem_stock,
  aes(x = stock_kg_m2, y = ecosystem, fill = dataset)
) +
  geom_boxplot(outlier.alpha = 0.15) +
  geom_text(
    data = ecosystem_counts,
    aes(x = Inf, y = ecosystem, label = paste0("n = ", n)),
    inherit.aes = FALSE, hjust = 1.05, size = 3
  ) +
  facet_wrap(~ interval, scales = "free_x", ncol = 2) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  labs(
    title = "Reference carbon stocks by ecosystem",
    x = "Carbon stock (kg C/m²)",
    y = NULL,
    fill = "Reference collection"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

ecosystem_plot <- ecosystem_plot +
  geom_vline(
    data = field,
    aes(xintercept = value, colour = core_id),
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  labs(colour = "Community core")

print(ecosystem_plot)

ggsave(
  "outputs/script1/carbon_stock_by_ecosystem_horizontal.png",
  ecosystem_plot,
  width = 14, height = 9, dpi = 200
)


#====
#----- Summary statistics for Part 3 --

# Where does our data fit in with larger soil and sediment data?


#====
#----- Summary statistics for Part 3 --
# Where does our data fit within larger soil and sediment datasets?

depth_order <- c("0-30", "0-50", "0-100", "total")

field_for_summary <- field %>%
  filter(is.finite(value)) %>%
  mutate(interval = as.character(interval))

reference_for_summary <- reference %>%
  filter(
    variable == "stock_kgm2",
    datum == "surface",
    unit == "kg C/m2",
    is.finite(value)
  ) %>%
  mutate(interval = as.character(interval))

if (!nrow(field_for_summary)) {
  stop("No valid field stocks are available to summarize.")
}

if (anyNA(field_for_summary$core_id) ||
    anyNA(field_for_summary$interval) ||
    anyDuplicated(field_for_summary[c("core_id", "interval")])) {
  stop("Each field core needs one unique record per depth group.")
}

reference_percentile <- function(stock, depth) {
  values <- reference_for_summary$value[
    reference_for_summary$interval == depth
  ]
  
  if (!length(values)) return(NA_real_)
  100 * mean(values <= stock)
}

part3_summary <- field_for_summary %>%
  group_by(interval) %>%
  summarise(
    field_cores = n(),
    mean_stock = mean(value),
    median_stock = median(value),
    sd_stock = if (n() > 1) sd(value) else NA_real_,
    minimum_stock = min(value),
    maximum_stock = max(value),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    reference_profiles = sum(
      reference_for_summary$interval == interval
    ),
    median_percentile = reference_percentile(median_stock, interval),
    minimum_percentile = reference_percentile(minimum_stock, interval),
    maximum_percentile = reference_percentile(maximum_stock, interval)
  ) %>%
  ungroup() %>%
  arrange(match(interval, depth_order))

part3_chart <- ggplot(part3_summary, aes(y = interval)) +
  annotate(
    "rect",
    xmin = 25, xmax = 75, ymin = -Inf, ymax = Inf,
    fill = "grey93"
  ) +
  geom_vline(
    xintercept = 50, linetype = "dashed", colour = "grey60"
  ) +
  geom_segment(
    aes(
      x = minimum_percentile,
      xend = maximum_percentile,
      yend = interval
    ),
    linewidth = 1.5, colour = "#34778c", na.rm = TRUE
  ) +
  geom_point(
    aes(x = median_percentile),
    shape = 18, size = 4, colour = "#163d4a", na.rm = TRUE
  ) +
  geom_text(
    data = filter(part3_summary, reference_profiles == 0),
    aes(x = 50, label = "No matching reference"),
    colour = "grey40"
  ) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_y_discrete(
    limits = rev(depth_order[depth_order %in% part3_summary$interval])
  ) +
  labs(
    title = "Where our field stocks sit within the reference data",
    x = "Percentile among available soil and sediment reference records",
    y = NULL,
    caption = paste(
      "Diamond = field median; line = field-core range; grey = reference middle 50%.",
      "Reporting groups retain the agreed depth caps; sampled depths can differ.",
      sep = "\n"
    )
  ) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank())

print(part3_chart)
print(part3_summary)

# Prefer 0–30 cm for the reporting sentence.
# Otherwise use total, then another available group.
summary_depth <- if ("0-30" %in% part3_summary$interval) {
  "0-30"
} else if ("total" %in% part3_summary$interval) {
  "total"
} else {
  part3_summary$interval[1]
}

s <- filter(part3_summary, interval == summary_depth)

depth_label <- if (summary_depth == "total") {
  "total sampled-profile"
} else {
  paste0(summary_depth, " cm")
}

part3_sentence <- sprintf(
  paste0(
    "For %s, field stocks (n = %d) averaged %.2f kg C/m² ",
    "(median %.2f; range %.2f–%.2f)"
  ),
  depth_label, s$field_cores, s$mean_stock,
  s$median_stock, s$minimum_stock, s$maximum_stock
)

if (s$reference_profiles > 0) {
  part3_sentence <- paste0(part3_sentence, sprintf(
    paste0(
      "; %.0f%% of %s available soil and sediment reference records ",
      "in this reporting group had stocks at or below that median."
    ),
    s$median_percentile,
    format(s$reference_profiles, big.mark = ",", trim = TRUE)
  ))
} else {
  part3_sentence <- paste0(
    part3_sentence, "; no matching reference records were available."
  )
}

cat(part3_sentence, "\n")

dir.create("outputs/script1", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "outputs/script1/part3_reference_position.png",
  part3_chart, width = 8.5, height = 3.5, dpi = 200
)

readr::write_csv(
  part3_summary, "outputs/script1/part3_summary_statistics.csv"
)

writeLines(
  part3_sentence, "outputs/script1/part3_summary_sentence.txt"
)

part3_summary
















