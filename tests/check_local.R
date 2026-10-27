source("R/local.R")
# Synthetic fixtures test real-data interfaces; none become workshop evidence.
out <- tempfile("interface-test-",tmpdir="outputs");dir.create(out)
inputs <- tempfile();dir.create(inputs)
settings <- list(field_path="data/input/field_data_template.csv",aoi_path="data/aoi/study_area.geojson",
  output_dir=out,cellsize_m=1000)
r <- run_local_workshop(settings)
stopifnot(all(read.csv(file.path(out,"local_stage_status.csv"))$status[-1]=="Not run: input path is NULL"))
z <- read.csv(file.path(out,"local_grid_template.csv"))
# Choose arbitrary ecosystem labels only for this software test.
z$ecosystem <- r$maps$cores$ecosystem[1]
z$wetness <- sin(z$x_m/3000);z$elevation_m <- z$y_m/1000
z$prior_mean_kg_m2 <- 15;z$prior_sd_kg_m2 <- 5;z$new_sd_kg_m2 <- 3
path <- file.path(inputs,"cells.csv");write.csv(z,path,row.names=FALSE)
points <- .grid_points(r$maps$cells)
keep <- lengths(sf::st_intersects(points,r$maps$aoi))>0
idx <- which(keep)[seq(1,sum(keep),length.out=60)]
obs <- z[idx, ]; xy <- sf::st_coordinates(sf::st_transform(points[idx,],4326))
obs$core_id <- paste0("TEST-",seq_len(nrow(obs)));obs$longitude <- xy[,1];obs$latitude <- xy[,2]
obs$stock_0_30_kg_m2 <- 20+2*obs$wetness-.1*obs$elevation_m
cut <- quantile(obs$x_m,.7)
obs$split <- ifelse(obs$x_m>cut,"Spatial holdout","Training")
obspath <- file.path(inputs,"obs.csv");write.csv(obs,obspath,row.names=FALSE)
settings$ecosystem_cells <- path;settings$predictor_cells <- path;settings$model_observations <- obspath
settings$prior_cells <- path;settings$independent_errors_confirmed <- TRUE;settings$holdout_buffer_m <- 10
r <- run_local_workshop(settings)
stopifnot(all(grepl("Complete",read.csv(file.path(out,"local_stage_status.csv"))$status)),
  read.csv(file.path(out,"local_model_metrics.csv"))$RMSE<1e-8)
totals <- read.csv(file.path(out,"local_area_totals.csv"))
stopifnot(abs(sum(totals$area_m2)/as.numeric(sf::st_area(r$maps$aoi))-1)<1e-7)
# Changed AOI geometry cannot masquerade as an old grid through recycled cell IDs.
z$x_m[1] <- z$x_m[1]+100
write.csv(z,path,row.names=FALSE)
stopifnot(inherits(try(.read_cell_table(path,r$maps$cells,"ecosystem"),silent=TRUE),"try-error"))
settings$field_path <- file.path(inputs,"missing.csv")
stopifnot(inherits(try(run_course(settings$field_path,settings$aoi_path,output_dir=out),silent=TRUE),"try-error"),
 !file.exists(file.path(out,"local_bayesian.csv")),!file.exists(file.path(out,"maps/idw_stock_map.png")),
 !file.exists(file.path(out,"part1/comparison/stock_by_dataset.png")))
unlink(out,recursive=TRUE);unlink(inputs,recursive=TRUE)
message("Local interfaces passed: explicit skips, all optional inputs, geometry alignment, held-out regression, Bayesian outputs and failed-rerun cleanup.")
