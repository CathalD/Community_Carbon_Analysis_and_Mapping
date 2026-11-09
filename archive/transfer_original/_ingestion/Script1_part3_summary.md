# Part 3 summary statistics and reference position

Paste this after the existing Part 3 comparisons and before the export section. It uses the current field and reference tables created earlier in Script 1. Re-running Script 1 on a new workbook in the same input format rebuilds the summary. No example core IDs, sample counts, or stock values are hard-coded.

Outputs: per-group mean, median, SD, range, field counts, reference counts, median/range percentiles; a small percentile-position chart; and one reporting sentence. SD is unavailable for a single observation. The sentence prefers 0–30 cm; if absent, it uses total sampled stock or another available group. Missing stocks are omitted, duplicate core/group entries are rejected, absent reference matches remain unavailable, and empty field data produces a clear error.

Reference percentiles pool the available profile records in the same reporting group. Collections are not equally weighted and records are not asserted to be independent sites. This is a broad reference comparison, not an ecosystem-matched assessment or a measure of habitat condition. The previously agreed capped-depth groups are retained, including partial sampled depths. The line in the small chart shows the field range, not a confidence interval.

```r
#====
#----- Summary statistics for Part 3 --
# Where does our data fit within larger soil and sediment datasets?

depth_order <- c("0-30", "0-50", "0-100", "total")
# Use the current observations, not previously exported example results.
field_for_summary <- field %>%
  filter(is.finite(value)) %>%
  mutate(interval = as.character(interval))
reference_for_summary <- reference %>%
  filter(variable == "stock_kgm2", datum == "surface",
         unit == "kg C/m2", is.finite(value)) %>%
  mutate(interval = as.character(interval))

if (!nrow(field_for_summary)) stop("No valid field stocks are available to summarize.")
if (anyNA(field_for_summary$core_id) || anyNA(field_for_summary$interval) ||
    anyDuplicated(field_for_summary[c("core_id", "interval")])) {
  stop("Each field core needs one unique record per depth group.")
}

reference_percentile <- function(stock, depth) {
  values <- reference_for_summary$value[reference_for_summary$interval == depth]
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
    reference_profiles = sum(reference_for_summary$interval == interval),
    median_percentile = reference_percentile(median_stock, interval),
    minimum_percentile = reference_percentile(minimum_stock, interval),
    maximum_percentile = reference_percentile(maximum_stock, interval)
  ) %>%
  ungroup() %>%
  arrange(match(interval, depth_order))

part3_chart <- ggplot(part3_summary, aes(y = interval)) +
  annotate("rect", xmin = 25, xmax = 75, ymin = -Inf, ymax = Inf,
           fill = "grey93") +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey60") +
  geom_segment(aes(x = minimum_percentile, xend = maximum_percentile,
                   yend = interval), linewidth = 1.5, colour = "#34778c", na.rm = TRUE) +
  geom_point(aes(x = median_percentile), shape = 18, size = 4,
             colour = "#163d4a", na.rm = TRUE) +
  geom_text(data = filter(part3_summary, reference_profiles == 0),
            aes(x = 50, label = "No matching reference"), colour = "grey40") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_discrete(limits = rev(depth_order[depth_order %in% part3_summary$interval])) +
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

# Prefer 0–30 cm; otherwise use total, then the first available group.
summary_depth <- if ("0-30" %in% part3_summary$interval) "0-30" else
  if ("total" %in% part3_summary$interval) "total" else part3_summary$interval[1]
s <- filter(part3_summary, interval == summary_depth)
depth_label <- if (summary_depth == "total") "total sampled-profile" else
  paste0(summary_depth, " cm")
part3_sentence <- sprintf(
  "For %s, field stocks (n = %d) averaged %.2f kg C/m² (median %.2f; range %.2f–%.2f)",
  depth_label, s$field_cores, s$mean_stock, s$median_stock,
  s$minimum_stock, s$maximum_stock
)
if (s$reference_profiles > 0) {
  part3_sentence <- paste0(part3_sentence, sprintf(
    "; %.0f%% of %s available soil and sediment reference records in this reporting group had stocks at or below that median.",
    s$median_percentile, format(s$reference_profiles, big.mark = ",", trim = TRUE)
  ))
} else {
  part3_sentence <- paste0(part3_sentence, "; no matching reference records were available.")
}
cat(part3_sentence, "\n")

dir.create("outputs/script1", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/script1/part3_reference_position.png", part3_chart,
       width = 8.5, height = 3.5, dpi = 200)
readr::write_csv(part3_summary, "outputs/script1/part3_summary_statistics.csv")
writeLines(part3_sentence, "outputs/script1/part3_summary_sentence.txt")
```

## One-line report template

For [depth group], field stocks (n = [valid cores]) averaged [mean] kg C/m² (median [median]; range [minimum]–[maximum]); [percent]% of [reference count] available soil and sediment reference records in this reporting group had stocks at or below that median.

The saved output is a preview using the current six simulated cores. Scripts 1–5 were not edited.
