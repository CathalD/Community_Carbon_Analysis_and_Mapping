
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
workbook_path = "data/input/community_carbon_outputs_updated.xlsm"

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
CuOrgC_plot

# And some summary stats to go with it:
# Calculate summary stats per Core ID
library(dplyr)

# Calculate summary stats per Core ID with depth-specific carbon stocks
core_summary <- sample_sheet %>%
  group_by(`Core ID`) %>%
  summarise(
    `Total Sections`      = n(),
    `Max Depth (cm)`      = max(`Bottom depth (cm)`, na.rm = TRUE),
    `Avg Bulk Density`    = mean(`Bulk density (g/cm³)`, na.rm = TRUE),
    `Avg Carbon (%)`      = mean(`Organic carbon (%)`, na.rm = TRUE),
    
    # Carbon Stocks at specific depth horizons
    `C Stock 0-30cm`      = sum(`Stock, full thickness (kg C/m²)`[`Bottom depth (cm)` <= 30], na.rm = TRUE),
    `C Stock 0-50cm`      = sum(`Stock, full thickness (kg C/m²)`[`Bottom depth (cm)` <= 50], na.rm = TRUE),
    `C Stock 0-100cm`     = sum(`Stock, full thickness (kg C/m²)`[`Bottom depth (cm)` <= 100], na.rm = TRUE),
    
    `Total Carbon Stock`  = sum(`Stock, full thickness (kg C/m²)`, na.rm = TRUE)
  ) %>%
  # Round all numeric metrics to 2 decimal places
  mutate(across(where(is.numeric), ~ round(., 2)))

# View the new summary table with depth metrics
print(core_summary)


# ========================================================

# ---------- Step 3 ----------------

# ==================================

# Bring in the external data


# Plot the external data

# Overlay our data on top of it




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













