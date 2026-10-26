# Compare recalculated community values with precomputed SOD reference summaries.
# The download contains aggregate statistics, not the underlying reference profiles.
run_reference_comparison <- function(config = list(), cores = NULL, layers = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 before this step.")
  default <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x
  output_root <- default(config$output_dir, "outputs")
  out <- file.path(output_root, "part1", "comparison")
  reference_dir <- default(config$reference_dir, "data/reference")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cores)) cores <- read.csv(file.path(output_root, "part1", "core_stocks.csv"), stringsAsFactors = FALSE)
  if (is.null(layers)) layers <- read.csv(file.path(output_root, "part1", "layer_stocks.csv"), stringsAsFactors = FALSE)
  required <- c("core_id", "top_cm", "bottom_cm", "bd_g_cm3", "carbon_used_pct")
  if (!all(required %in% names(layers))) stop("Layer outputs are missing: ", paste(setdiff(required, names(layers)), collapse = ", "))
  if (!"core_id" %in% names(cores) || anyDuplicated(cores$core_id)) stop("Core outputs require unique core_id values.")
  if (any(!layers$core_id %in% cores$core_id)) stop("Some layers have no matching core record.")
  if (any(!is.finite(layers$top_cm) | !is.finite(layers$bottom_cm) | layers$top_cm < 0 | layers$bottom_cm <= layers$top_cm)) stop("Check layer depths before comparing.")
  for (z in split(layers, layers$core_id)) {
    z <- z[order(z$top_cm), ]
    if (nrow(z) > 1 && any(z$top_cm[-1] < head(z$bottom_cm, -1) - 1e-8)) stop("Overlapping layers in ", z$core_id[1])
  }
  dimensions <- c("dataset", "study_group", "ecosystem", "general_region", "source_region")
  summary_columns <- c("dataset", "variable", "interval", "datum", "unit", "group", "n_profile_records", "recorded_coordinate_pairs", "mean", "median", "sample_sd", "cv_pct", "q25", "q75", "min", "max")
  summaries <- setNames(lapply(dimensions, function(dimension) {
    path <- file.path(reference_dir, paste0("summary_by_", dimension, ".csv.gz"))
    if (!file.exists(path)) stop("Reference summary missing: ", path)
    con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
    tab <- read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(summary_columns %in% names(tab))) stop("Unexpected summary columns in ", path)
    if (anyDuplicated(tab[c("dataset", "variable", "interval", "datum", "group")])) stop("Duplicate reference summary groups in ", path)
    if (any(tab$n_profile_records < 1 | !is.finite(tab$mean) | !is.finite(tab$median))) stop("Invalid reference summaries in ", path)
    tab$group[is.na(tab$group) | !nzchar(tab$group)] <- "Not recorded"
    tab
  }), dimensions)
  # Recalculate each community quantity from its imported layers. A slice that
  # crosses a reporting boundary is allocated by overlap; no gaps are filled.
  windows <- unique(summaries$dataset[c("variable", "interval", "datum", "unit")])
  datum <- default(config$depth_datum, "surface")
  if (!datum %in% c("surface", "mineral")) stop("config$depth_datum must be surface or mineral.")
  windows <- windows[windows$datum == datum, ]
  origin <- as.character(default(config$data_origin, "Community input; origin not specified"))[1]
  community <- do.call(rbind, lapply(split(layers, layers$core_id), function(z) {
    z <- z[order(z$top_cm), ]
    core <- cores[match(z$core_id[1], cores$core_id), ]
    rows <- lapply(seq_len(nrow(windows)), function(i) {
      w <- windows[i, ]
      bounds <- as.numeric(strsplit(w$interval, "-", fixed = TRUE)[[1]])
      overlap <- pmax(0, pmin(z$bottom_cm, bounds[2]) - pmax(z$top_cm, bounds[1]))
      keep <- overlap > 0
      clipped_top <- pmax(z$top_cm[keep], bounds[1])
      clipped_bottom <- pmin(z$bottom_cm[keep], bounds[2])
      complete <- any(keep) && abs(sum(overlap) - diff(bounds)) < 1e-8 &&
        abs(min(clipped_top) - bounds[1]) < 1e-8 && abs(max(clipped_bottom) - bounds[2]) < 1e-8 &&
        (length(clipped_top) == 1 || all(abs(clipped_top[-1] - head(clipped_bottom, -1)) < 1e-8))
      property <- switch(w$variable,
        soc_gkg = z$carbon_used_pct * 10,
        bd = z$bd_g_cm3,
        stock_kgm2 = z$bd_g_cm3 * z$carbon_used_pct * 0.1,
        stop("Unknown comparison variable: ", w$variable))
      valid <- complete && all(is.finite(property[keep]))
      value <- if (!valid) NA_real_ else if (w$variable == "stock_kgm2") sum(property[keep] * overlap[keep]) else sum(property[keep] * overlap[keep]) / diff(bounds)
      ecosystem <- if ("ecosystem" %in% names(core)) as.character(core$ecosystem[1]) else "Not recorded"
      group <- if ("group" %in% names(core)) as.character(core$group[1]) else "Not recorded"
      if (is.na(ecosystem) || !nzchar(ecosystem)) ecosystem <- "Not recorded"
      if (is.na(group) || !nzchar(group)) group <- "Not recorded"
      data.frame(core_id = core$core_id, ecosystem = ecosystem, community_group = group,
        data_origin = origin, variable = w$variable, interval = w$interval, datum = w$datum,
        unit = w$unit, value = value, complete_interval = complete,
        eligible = valid, note = if (valid) "" else if (!complete) "Incomplete sampled interval; no extrapolation" else "Missing measurement", stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }))
  rownames(community) <- NULL
  write.csv(community, file.path(out, "community_matched_depth_values.csv"), row.names = FALSE, na = "")
  all_summaries <- do.call(rbind, lapply(dimensions, function(dimension) {
    z <- summaries[[dimension]]; z$dimension <- dimension; z
  }))
  # Join on quantity, units, interval and datum only. Ecosystem/geographic
  # similarity is not asserted: these are broad descriptive comparisons.
  comparisons <- merge(community[community$eligible, ], all_summaries,
    by = c("variable", "interval", "datum", "unit"), all = FALSE, sort = FALSE)
  comparisons$difference_from_reference_mean <- comparisons$value - comparisons$mean
  comparisons$difference_from_reference_median <- comparisons$value - comparisons$median
  comparisons$pct_difference_from_reference_median <- ifelse(comparisons$median != 0,
    100 * (comparisons$value - comparisons$median) / comparisons$median, NA_real_)
  comparisons$within_reference_iqr <- comparisons$value >= comparisons$q25 & comparisons$value <= comparisons$q75
  comparisons$comparison_scope <- "Same quantity, units, depth and datum; broad context, not a matched ecosystem test"
  write.csv(comparisons, file.path(out, "community_vs_reference_groups.csv"), row.names = FALSE, na = "")
  availability <- summaries$dataset[c("dataset", "variable", "interval", "datum", "unit", "n_profile_records", "recorded_coordinate_pairs")]
  write.csv(availability, file.path(out, "reference_availability.csv"), row.names = FALSE)
  # Export the complete aggregate tables as well as the smaller figure views.
  for (dimension in dimensions) write.csv(summaries[[dimension]],
    file.path(out, paste0("reference_summary_by_", dimension, ".csv")), row.names = FALSE, na = "")
  figure_data <- list()
  plot_summary <- function(dimension, variable, interval, filename, title, x_label, top_n = Inf) {
    tab <- summaries[[dimension]]
    tab <- tab[tab$variable == variable & tab$interval == interval & tab$datum == datum, ]
    tab <- tab[order(-tab$n_profile_records, tab$dataset, tab$group), ]
    if (nrow(tab) > top_n) tab <- head(tab, top_n)
    if (!nrow(tab)) return(invisible(NULL))
    label <- if (dimension == "dataset") tab$dataset else if (dimension == "study_group") tab$group else paste0(tab$group, " [", tab$dataset, "]")
    tab$display_label <- label
    tab$display_label <- factor(tab$display_label, levels = rev(unique(label)))
    obs <- community[community$eligible & community$variable == variable & community$interval == interval & community$datum == datum, ]
    p <- ggplot2::ggplot(tab, ggplot2::aes(y = display_label)) +
      ggplot2::geom_segment(ggplot2::aes(x = q25, xend = q75, yend = display_label), linewidth = 2.2, colour = "#9EB7C5") +
      ggplot2::geom_point(ggplot2::aes(x = median, shape = "Reference median"), size = 2.7, colour = "#203B4B") +
      ggplot2::geom_point(ggplot2::aes(x = mean, shape = "Reference mean"), size = 2.9, colour = "#203B4B", fill = "white") +
      ggplot2::scale_shape_manual(values = c("Reference median" = 16, "Reference mean" = 23), name = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom",
        axis.text.y = ggplot2::element_text(size = 9)) +
      ggplot2::labs(x = x_label, y = NULL)
    if (nrow(obs)) p <- p + ggplot2::geom_vline(data = obs,
      ggplot2::aes(xintercept = value, colour = core_id), linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
      ggplot2::labs(colour = "Community core") + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))
    ggplot2::ggsave(file.path(out, filename), p, width = 12,
      height = max(5, min(13, 2.8 + 0.34 * nrow(tab))), dpi = 160, limitsize = FALSE)
    tab$figure <- filename; figure_data[[filename]] <<- tab
    invisible(p)
  }
  plot_summary("dataset", "stock_kgm2", "0-30", "stock_by_dataset.png", "Carbon stock by reference collection", "Soil organic carbon, 0-30 cm (kg C/m2)")
  plot_summary("ecosystem", "stock_kgm2", "0-30", "stock_by_ecosystem.png", "Carbon stock by recorded ecosystem", "Soil organic carbon, 0-30 cm (kg C/m2)")
  plot_summary("general_region", "stock_kgm2", "0-30", "stock_by_region.png", "Carbon stock by general region", "Soil organic carbon, 0-30 cm (kg C/m2)")
  plot_summary("study_group", "stock_kgm2", "0-30", "stock_by_study.png", "Carbon stock by study or contributing collection", "Soil organic carbon, 0-30 cm (kg C/m2)", top_n = 20)
  plot_summary("source_region", "stock_kgm2", "0-30", "stock_by_source_region.png", "Carbon stock by province or source region", "Soil organic carbon, 0-30 cm (kg C/m2)", top_n = 20)
  plot_summary("dataset", "soc_gkg", "0-5", "concentration_by_dataset.png", "Organic carbon concentration by reference collection", "Organic carbon, 0-5 cm (g C/kg)")
  plot_summary("dataset", "bd", "0-5", "density_by_dataset.png", "Dry bulk density by reference collection", "Dry bulk density, 0-5 cm (g/cm3)")
  displayed <- do.call(rbind, figure_data); rownames(displayed) <- NULL
  write.csv(displayed, file.path(out, "figure_reference_groups.csv"), row.names = FALSE, na = "")
  source_files <- file.path(reference_dir, paste0("summary_by_", dimensions, ".csv.gz"))
  write.csv(data.frame(file = basename(source_files), md5 = unname(tools::md5sum(source_files))),
    file.path(out, "reference_input_checksums.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))
  message("Reference comparisons saved to ", out,
    ". Community values were recalculated; reference statistics use the supplied aggregate snapshot.")
  invisible(list(community_values = community, summaries = summaries, comparisons = comparisons,
    availability = availability, figure_groups = displayed, output_dir = out))
}
