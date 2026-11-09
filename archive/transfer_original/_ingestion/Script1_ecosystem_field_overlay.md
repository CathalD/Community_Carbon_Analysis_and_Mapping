# Add the field-core overlay

Paste this immediately before `print(ecosystem_plot)` and the existing `ggsave()` call. It uses the `field` table already created in Part 3 and adds only the vertical dashed lines and their legend. Reference data, depth groups and styling are retained.

```r
ecosystem_plot <- ecosystem_plot +
  geom_vline(
    data = field,
    aes(xintercept = value, colour = core_id),
    linetype = "dashed", linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  labs(colour = "Community core")
```

The saved chart includes 6 field lines at 0–30, 5 at 0–50, 4 at 0–100, and 6 at total. Field lines are broad reference overlays across ecosystem rows, with no habitat assignment implied. Scripts 1–5 and the reference CSV were not edited.
