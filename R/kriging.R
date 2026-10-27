# Ordinary kriging demonstration on a separate stationary Gaussian simulation.
# Nothing in this file estimates covariance from the four workbook cores.
source("R/spatial.R")
.fit_teaching_variogram <- function(points) {
  xy <- sf::st_coordinates(points)
  span <- sqrt(sum((apply(xy,2,max)-apply(xy,2,min))^2))
  empirical <- gstat::variogram(stock_0_30_kg_m2 ~ 1, points,
    cutoff=span*.6, width=span*.6/12)
  empirical <- empirical[empirical$np >= 10, ]
  if (nrow(empirical)<5) stop("Too few supported lag bins to fit a teaching variogram.")
  v <- stats::var(points$stock_0_30_kg_m2)
  # Restrict the search to explicit plausible ranges; no unbounded three-parameter fit.
  candidates <- lapply(span*c(.06,.10,.16,.25,.4), function(r) {
    suppressWarnings(gstat::fit.variogram(empirical,
      gstat::vgm(v*.8,"Exp",r,v*.2), fit.ranges=FALSE, fit.method=7))
  })
  valid <- vapply(candidates,function(m) !isTRUE(attr(m,"singular")) &&
    all(is.finite(m$psill)) && all(m$psill>=0) && m$psill[m$model=="Exp"]>0 &&
    is.finite(attr(m,"SSErr")),logical(1))
  if (!any(valid)) stop("Variogram fit failed. Inspect spatial coverage and lag bins; no fallback map was made.")
  candidates <- candidates[valid]
  m <- candidates[[which.min(vapply(candidates,function(z) attr(z,"SSErr"),numeric(1)))]]
  list(model=m, empirical=empirical)
}
run_kriging_demo <- function(config, area_path="outputs/part2/area/area_analysis.rds",
  out="outputs/part2/kriging") {
  if (!isTRUE(config$illustrative)) stop("Kriging demonstration requires illustrative = TRUE; real covariance needs a separate assessment.")
  a <- readRDS(area_path)
  set.seed(config$seed + 1L)
  locations <- sf::st_sf(geometry=sf::st_sample(a$aoi, size=240, exact=TRUE))
  xy <- sf::st_coordinates(locations)
  span <- sqrt(sum((apply(xy,2,max)-apply(xy,2,min))^2))
  truth <- gstat::vgm(psill=9, model="Exp", range=span/6, nugget=1)
  generator <- gstat::gstat(formula=z~1, locations=~x+y, dummy=TRUE,
    beta=18, model=truth, nmax=240)
  sim <- predict(generator, newdata=data.frame(x=xy[,1],y=xy[,2]), nsim=1, debug.level=0)
  locations$stock_0_30_kg_m2 <- sim$sim1
  locations$record_id <- seq_len(nrow(locations))
  locations$x_m <- xy[,1]; locations$y_m <- xy[,2]
  cutoff <- unname(stats::quantile(xy[,1],.75))
  buffer <- diff(range(xy[,1]))*.05
  locations$split <- ifelse(xy[,1]>cutoff,"Spatial holdout",
    ifelse(xy[,1]<cutoff-buffer,"Training","Buffer"))
  locations$data_origin <- "SIMULATED stationary Gaussian stock; not field observations"
  train <- locations[locations$split=="Training", ]; test <- locations[locations$split=="Spatial holdout", ]
  fitted <- .fit_teaching_variogram(train)
  prediction <- gstat::krige(stock_0_30_kg_m2~1,train,test,model=fitted$model,debug.level=0)
  test$predicted <- prediction$var1.pred; test$variance <- prediction$var1.var
  test$idw_predicted <- .gstat_idw(test,train)
  test$residual <- test$stock_0_30_kg_m2-test$predicted
  if(any(!is.finite(test$variance) | test$variance<=0)) stop("Invalid held-out kriging variance.")
  metrics <- data.frame(method=c("Ordinary kriging","IDW power 2"),
    n_training=nrow(train),n_holdout=nrow(test), buffer_m=buffer,
    RMSE=c(sqrt(mean(test$residual^2)),sqrt(mean((test$stock_0_30_kg_m2-test$idw_predicted)^2))),
    MAE=c(mean(abs(test$residual)),mean(abs(test$stock_0_30_kg_m2-test$idw_predicted))),
    mean_standardised_error=c(mean(test$residual/sqrt(test$variance)),NA),
    coverage95=c(mean(abs(test$residual)<=stats::qnorm(.975)*sqrt(test$variance)),NA),
    scope="One buffered eastern simulation holdout; conditional model variance; not field validation")
  final <- .fit_teaching_variogram(locations)
  predicted <- gstat::krige(stock_0_30_kg_m2~1,locations,.grid_points(a$cells),
    model=final$model,debug.level=0)
  cells <- a$cells
  cells$kriging_mean_kg_m2 <- predicted$var1.pred
  cells$kriging_sd_kg_m2 <- sqrt(pmax(0,predicted$var1.var))
  if(any(!is.finite(cells$kriging_mean_kg_m2)) || any(predicted$var1.var < -1e-8)) stop("Invalid kriging prediction.")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(nm in c("simulation","holdout","metrics","training_variogram","fitted_variogram","generating_variogram","grid_values")) {
    z <- switch(nm, simulation=sf::st_drop_geometry(locations), holdout=sf::st_drop_geometry(test),
      metrics=metrics, training_variogram=fitted$empirical,
      fitted_variogram=rbind(cbind(stage="Training only",fitted$model),cbind(stage="Final all simulation",final$model)),
      generating_variogram=truth, grid_values=sf::st_drop_geometry(cells))
    write.csv(z,file.path(out,paste0(nm,".csv")),row.names=FALSE)
  }
  lines <- gstat::variogramLine(fitted$model,maxdist=max(fitted$empirical$dist))
  p <- ggplot2::ggplot(fitted$empirical,ggplot2::aes(dist,gamma)) +
    ggplot2::geom_point(ggplot2::aes(size=np),colour="#317873") +
    ggplot2::geom_line(data=lines) + ggplot2::theme_bw() +
    ggplot2::labs(x="Separation (m)",y="Semivariance (kg C/m²)²",size="Pairs")
  ggplot2::ggsave(file.path(out,"variogram.png"),p,width=7,height=5,dpi=180)
  .save_map_panels(list(.stock_map(cells,a$aoi,"kriging_mean_kg_m2",panel="Simulated kriging mean"),
    .stock_map(cells,a$aoi,"kriging_sd_kg_m2",panel="Conditional kriging SD")),
    file.path(out,"kriging_maps.png"),width=14,height=6)
  for(nm in c("kriging_mean_kg_m2","kriging_sd_kg_m2"))
    .write_raster(cells,nm,file.path(out,paste0(nm,".tif")),a$cellsize_m)
  saveRDS(list(cells=cells,training_model=fitted$model,final_model=final$model,
    simulation=locations,holdout=test,metrics=metrics),file.path(out,"kriging_analysis.rds"))
  .spatial_metadata(out,"Ordinary point kriging on a separate simulation",
    "prediction=sum(lambda_i*z_i); sum(lambda_i)=1; gamma(h)=nugget+partial_sill*(1-exp(-h/range)) for h>0",
    "gstat::gstat; gstat::variogram; gstat::fit.variogram; gstat::krige; sf; terra; ggplot2",
    paste("240 simulated point observations; stationary isotropic exponential covariance, constant unknown mean.",
      "Range selected from five candidates using training variogram SSE; training-only fit for buffered holdout.",
      "Final variogram refitted after assessment. Gaussian point support, including nugget variation.",
      "Conditional variance excludes fitted-parameter uncertainty. Pixel values are point predictions, not block means.",
      "Four workbook cores are not used. No AOI stock total or real-map accuracy is claimed."))
  invisible(metrics)
}
