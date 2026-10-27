# Optional real-data interfaces. No missing predictors, ecosystems or SDs are invented.
source("R/course.R")
source("R/compare.R")
source("R/spatial.R")
.read_cell_table <- function(path, cells, required) {
  if (!file.exists(path)) stop("Missing local input: ", path)
  z <- read.csv(path, stringsAsFactors=FALSE)
  if (!all(c("cell_id",required) %in% names(z))) stop("Required columns: ",paste(c("cell_id",required),collapse=", "))
  if(anyDuplicated(z$cell_id) || !setequal(z$cell_id,cells$cell_id))
    stop("cell_id must match every current grid cell exactly once. Regenerate inputs after an AOI/resolution change.")
  z <- z[match(cells$cell_id,z$cell_id), ,drop=FALSE]
  # Exact coordinates and clipped areas guard against reusing IDs from an older AOI.
  for(nm in c("x_m","y_m","area_m2")) {
    if(!nm %in% names(z) || any(!is.finite(z[[nm]])) || any(abs(z[[nm]]-cells[[nm]])>1e-5))
      stop("Grid geometry mismatch in ",nm,". Use the current local_grid_template.csv.")
  }
  z
}
run_local_workshop <- function(settings) {
  out <- settings$output_dir
  if(is.null(out) || !grepl("^outputs/",out) || grepl("\\.\\.",out)) stop("Use an output_dir below outputs/.")
  # Reserve only these named advanced outputs; preserve all learner-added files.
  owned <- c("local_grid_template.csv","local_grid.gpkg","local_stage_status.csv",
    "local_area_totals.csv","local_group_map.png","local_model_map.png","local_holdout.csv",
    "local_model_metrics.csv","local_model_coefficients.csv","local_bayesian.csv","local_bayesian_map.png")
  unlink(file.path(out,owned))
  source("R/pipeline.R",local=TRUE)
  unlink(file.path(out,"part1/comparison",basename(.stage_files$comparison)))
  complete <- FALSE
  on.exit(if(!complete) unlink(file.path(out,owned)),add=TRUE)
  result <- run_course(settings$field_path,settings$aoi_path,output_dir=out,
    grid_cellsize_m=settings$cellsize_m)
  run_reference_comparison(list(output_dir=out,depth_datum="surface",data_origin="User-supplied observations"))
  cells <- result$maps$cells
  write.csv(sf::st_drop_geometry(cells)[c("cell_id","x_m","y_m","area_m2")],
    file.path(out,"local_grid_template.csv"),row.names=FALSE)
  sf::st_write(cells,file.path(out,"local_grid.gpkg"),delete_dsn=TRUE,quiet=TRUE)
  status <- data.frame(stage=c("Field, IDW and reference comparison","Ecosystem totals","Landscape regression","Bayesian update"),
    status=c("Complete",rep("Not run: input path is NULL",3)))
  if(!is.null(settings$ecosystem_cells)) {
    z <- .read_cell_table(settings$ecosystem_cells,cells,"ecosystem")
    if(anyNA(z$ecosystem) || any(!nzchar(z$ecosystem))) stop("Every cell needs an ecosystem label from a real map.")
    core <- result$maps$cores
    means <- tapply(core$stock_0_30_kg_m2,core$ecosystem,mean)
    if(any(!z$ecosystem %in% names(means))) stop("A mapped ecosystem has no usable sampled core; do not fill it with another ecosystem mean.")
    cells$ecosystem <- z$ecosystem; cells$local_mean <- unname(means[z$ecosystem])
    totals <- do.call(rbind,lapply(split(sf::st_drop_geometry(cells),cells$ecosystem),function(v)
      data.frame(ecosystem=v$ecosystem[1],area_m2=sum(v$area_m2),mean_kg_m2=v$local_mean[1],
        n_cores=sum(core$ecosystem==v$ecosystem[1]),total_t_C=sum(v$local_mean*v$area_m2)/1000)))
    write.csv(totals,file.path(out,"local_area_totals.csv"),row.names=FALSE)
    ggplot2::ggsave(file.path(out,"local_group_map.png"),.stock_map(cells,result$maps$aoi,"local_mean"),width=7,height=6)
    status$status[2] <- "Complete; representation and ecosystem rasterisation remain user assumptions"
  }
  if(xor(is.null(settings$predictor_cells),is.null(settings$model_observations)))
    stop("Supply both predictor_cells and model_observations, or leave both NULL.")
  if(!is.null(settings$predictor_cells)) {
    predictors <- c("wetness","elevation_m")
    z <- .read_cell_table(settings$predictor_cells,cells,predictors)
    obs <- read.csv(settings$model_observations,stringsAsFactors=FALSE)
    required <- c("core_id","longitude","latitude","stock_0_30_kg_m2",predictors,"split")
    if(!all(required %in% names(obs))) stop("Model observations need: ",paste(required,collapse=", "))
    if(anyDuplicated(obs$core_id) || anyDuplicated(obs[c("longitude","latitude")])) stop("Model observations require distinct IDs and locations.")
    for(nm in c("longitude","latitude","stock_0_30_kg_m2",predictors))
      if(!is.numeric(obs[[nm]]) || any(!is.finite(obs[[nm]]))) stop("Missing/invalid observation ",nm)
    for(nm in predictors) if(!is.numeric(z[[nm]]) || any(!is.finite(z[[nm]]))) stop("Missing/invalid cell ",nm)
    if(any(abs(obs$longitude)>180 | abs(obs$latitude)>90) || any(obs$stock_0_30_kg_m2<0)) stop("Invalid model coordinates or negative stocks.")
    pts <- sf::st_transform(sf::st_as_sf(obs,coords=c("longitude","latitude"),crs=4326),sf::st_crs(cells))
    if(any(lengths(sf::st_intersects(pts,result$maps$aoi))==0)) stop("Model observations fall outside the AOI.")
    if(any(!obs$split %in% c("Training","Spatial holdout","Buffer"))) stop("Use Training, Spatial holdout or Buffer split labels.")
    train <- obs[obs$split=="Training", ]; test <- obs[obs$split=="Spatial holdout", ]
    if(nrow(train)<10 || nrow(test)<5) stop("Need at least ten training and five held-out observations; these are teaching safeguards, not adequacy guarantees.")
    if(is.null(settings$holdout_buffer_m) || !is.finite(settings$holdout_buffer_m) || settings$holdout_buffer_m<=0)
      stop("Set holdout_buffer_m to a scientifically justified positive separation.")
    separation <- .distance_matrix(pts[obs$split=="Training", ],pts[obs$split=="Spatial holdout", ])
    if(min(separation)<settings$holdout_buffer_m) stop("Training and holdout locations violate holdout_buffer_m.")
    fit <- stats::lm(stock_0_30_kg_m2~wetness+elevation_m,data=train)
    if(fit$rank<3) stop("Predictors are rank deficient in training observations.")
    test$predicted <- as.numeric(predict(fit,test))
    metrics <- data.frame(RMSE=sqrt(mean((test$stock_0_30_kg_m2-test$predicted)^2)),
      MAE=mean(abs(test$stock_0_30_kg_m2-test$predicted)),n_training=nrow(train),n_holdout=nrow(test),
      minimum_separation_m=min(separation))
    write.csv(test,file.path(out,"local_holdout.csv"),row.names=FALSE)
    write.csv(metrics,file.path(out,"local_model_metrics.csv"),row.names=FALSE)
    final <- stats::lm(stock_0_30_kg_m2~wetness+elevation_m,data=obs)
    cells$local_model <- as.numeric(predict(final,z))
    if(any(cells$local_model<0)) warning("Linear model predicts negative stock. Review model form; values were not truncated.")
    write.csv(data.frame(term=names(coef(final)),estimate=unname(coef(final))),file.path(out,"local_model_coefficients.csv"),row.names=FALSE)
    ggplot2::ggsave(file.path(out,"local_model_map.png"),.stock_map(cells,result$maps$aoi,"local_model"),width=7,height=6)
    status$status[3] <- "Complete; holdout assessed before final refit; predictor provenance must be reviewed"
  }
  if(!is.null(settings$prior_cells)) {
    required <- c("prior_mean_kg_m2","prior_sd_kg_m2","new_sd_kg_m2")
    z <- .read_cell_table(settings$prior_cells,cells,required)
    if(!isTRUE(settings$independent_errors_confirmed)) stop("Review prior/new error dependence and set independent_errors_confirmed=TRUE only if justified.")
    post <- normal_update(z$prior_mean_kg_m2,z$prior_sd_kg_m2,cells$idw_stock_0_30_kg_m2,z$new_sd_kg_m2)
    cells$local_posterior <- post$mean
    write.csv(cbind(z,new_mean_kg_m2=cells$idw_stock_0_30_kg_m2,post),file.path(out,"local_bayesian.csv"),row.names=FALSE)
    ggplot2::ggsave(file.path(out,"local_bayesian_map.png"),.stock_map(cells,result$maps$aoi,"local_posterior"),width=7,height=6)
    status$status[4] <- "Complete; supplied independent Gaussian errors; no whole-area uncertainty"
  }
  write.csv(status,file.path(out,"local_stage_status.csv"),row.names=FALSE)
  complete <- TRUE
  print(status,row.names=FALSE)
  invisible(result)
}
