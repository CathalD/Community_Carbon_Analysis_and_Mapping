# Preserve paired means/SDs. An isolated prior NoData cell must not silently drop
# a measured FS core. Only its conditioning prior can use a nearby valid pair;
# the original raster and its NoData mask remain unchanged.
fill_missing_core_priors <- function(native_prior, core_points, at_cores, ids,
                                    audit_file, max_distance_m = 500, allow_zero_mean = FALSE) {
  pts <- project(core_points, crs(native_prior))
  metric <- project(core_points, "EPSG:3978")
  audit <- data.frame(core_id = ids, method = "Original containing pixel",
    source_distance_m = 0, original_mean = at_cores$prior_mean, original_sd = at_cores$prior_sd,
    used_mean = at_cores$prior_mean, used_sd = at_cores$prior_sd)
  valid_mean <- function(x) is.finite(x) & if (allow_zero_mean) x >= 0 else x > 0
  bad <- which(!valid_mean(at_cores$prior_mean) |
    !is.finite(at_cores$prior_sd) | at_cores$prior_sd <= 0)
  for (i in bad) {
    window <- crop(native_prior, buffer(pts[i], width = max_distance_m * 2))
    pairs <- as.data.frame(window, xy = TRUE, na.rm = TRUE)
    pairs <- pairs[valid_mean(pairs$prior_mean) &
      is.finite(pairs$prior_sd) & pairs$prior_sd > 0, ]
    if (!nrow(pairs)) stop("No nearby paired prior for ", ids[i])
    coords <- crds(project(vect(pairs, geom = c("x", "y"), crs = crs(window)), crs(metric)))
    dist <- sqrt(rowSums(sweep(coords, 2, crds(metric)[i,], "-")^2))
    j <- which.min(dist)
    if (dist[j] > max_distance_m) stop("Paired prior gap exceeds ", max_distance_m, " m at ", ids[i])
    at_cores[i,] <- pairs[j, c("prior_mean", "prior_sd")]
    audit$method[i] <- "Nearest valid paired source-pixel centre within 500 m; conditioning only"
    audit$source_distance_m[i] <- dist[j]
    audit$used_mean[i] <- pairs$prior_mean[j]
    audit$used_sd[i] <- pairs$prior_sd[j]
  }
  write.csv(audit, audit_file, row.names = FALSE)
  at_cores
}
