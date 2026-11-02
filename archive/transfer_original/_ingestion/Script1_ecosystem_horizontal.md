# Horizontal ecosystem comparison

Paste this after the existing Part 3 comparison in Script 1. The existing four stock groups and measured-depth eligibility rules are retained. The CSV combines all 9,878 eligible stock records and includes ecosystem, source, source_url, original collection code, study group, actual sampled depths, and partial-depth flags. Records with ecosystem `Not recorded` remain in the CSV and are omitted from the ecosystem graphic. Riley has no computed stock values because its density units are unverified. The field example has unverified habitats, so field overlays are not assigned to ecosystem groups. No Script 1 source file was edited.

```r
# Add after the existing Part 3 comparison in Script 1.
ecosystem_stock <- readr::read_csv(
  "data/reference/ecosystem_carbon_stocks.csv",
  show_col_types = FALSE
) %>%
  filter(ecosystem != "Not recorded") %>%
  mutate(interval = factor(interval, levels = c("0-30", "0-50", "0-100", "total")))

ecosystem_counts <- ecosystem_stock %>%
  count(ecosystem, interval)

ecosystem_plot <- ggplot(
  ecosystem_stock,
  aes(x = stock_kg_m2, y = ecosystem, fill = dataset)
) +
  geom_boxplot(outlier.alpha = 0.15, linewidth = 0.35) +
  geom_text(
    data = ecosystem_counts,
    aes(x = Inf, y = ecosystem, label = paste0("n = ", n)),
    inherit.aes = FALSE, hjust = 1.05, size = 3
  ) +
  facet_wrap(~ interval, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("Janousek" = "#9ebdcc", "NRCan peat" = "#a4b998")) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  labs(
    title = "Reference carbon stocks by ecosystem",
    x = "Carbon stock (kg C/m²)", y = NULL, fill = "Reference collection",
    caption = paste(
      "Stocks are capped at the reporting depth; eligible shorter cores retain measured stock.",
      "Total covers each sampled profile. n = eligible profile records.",
      "Records without ecosystem labels remain in the CSV.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

print(ecosystem_plot)

dir.create("outputs/script1", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/script1/carbon_stock_by_ecosystem_horizontal.png",
       ecosystem_plot, width = 14, height = 9, dpi = 200)
```
