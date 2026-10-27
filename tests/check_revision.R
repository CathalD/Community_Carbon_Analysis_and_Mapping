source("R/spatial.R")
# Check gstat IDW against independent arithmetic, including exact observations.
p <- sf::st_as_sf(data.frame(x=c(0,100),y=c(0,0),stock_0_30_kg_m2=c(10,20)),coords=c("x","y"),crs=3978)
q <- sf::st_as_sf(data.frame(x=c(0,50,200),y=c(0,0,50)),coords=c("x","y"),crs=3978)
stopifnot(max(abs(.gstat_idw(q,p)-idw_values(.distance_matrix(q,p),p$stock_0_30_kg_m2)))<1e-10)
# Reject excessive pair allocation before calling the spatial library.
stopifnot(inherits(try(.distance_matrix(q[rep(1,2001),],p[rep(1,1000),]),silent=TRUE),"try-error"))
k <- readRDS("outputs/part2/kriging/kriging_analysis.rds")
tr <- k$simulation[k$simulation$split=="Training",]; te <- k$holdout
stopifnot(nrow(k$simulation)==240,all(tr$x_m < min(te$x_m)),
  !any(tr$record_id %in% te$record_id),all(is.finite(k$cells$kriging_sd_kg_m2)),
  all(k$cells$kriging_sd_kg_m2>=0),!isTRUE(attr(k$training_model,"singular")))
recheck <- gstat::krige(stock_0_30_kg_m2~1,tr,te,model=k$training_model,debug.level=0)
stopifnot(max(abs(recheck$var1.pred-te$predicted))<1e-10)
# Pipeline failure removes stale/partial products and downstream checkpoints, preserving notes.
source("R/pipeline.R")
old_files <- .stage_files
old_marker <- .stage_marker
sandbox <- tempfile();dir.create(sandbox)
.stage_marker <- file.path(sandbox,"run_complete.txt")
.stage_files <- lapply(.stage_files,function(x) file.path(sandbox,basename(x)))
for(f in unique(unlist(.stage_files))) writeLines("old",f)
writeLines("keep",file.path(sandbox,"notes.txt"))
try(run_stage("area", {writeLines("partial",.stage_files$area[1]);stop("test failure")}),silent=TRUE)
stopifnot(!any(file.exists(unlist(.stage_files[c("area","spatial","kriging","bayes")]))),
  file.exists(file.path(sandbox,"notes.txt")),all(file.exists(.stage_files$field)))
.stage_files <- old_files;.stage_marker <- old_marker;unlink(sandbox,recursive=TRUE)
# Existing portable cases plus an irregular polygon with a hole and a projected input.
# Avoid source(check_custom_area.R) here; it is already run in the main test suite.
source("config.R")
a <- .read_spatial_inputs(course_config)
b <- sf::st_bbox(a$aoi)
outer <- rbind(c(-10000,-10000),c(10000,-10000),c(7000,10000),c(-10000,4000),c(-10000,-10000))
hole <- rbind(c(-1000,-1000),c(-1000,1000),c(1000,1000),c(1000,-1000),c(-1000,-1000))
irregular <- sf::st_sf(geometry=sf::st_sfc(sf::st_polygon(list(outer,hole)),crs=sf::st_crs(a$aoi)))
g <- .area_grid(irregular,1000)
stopifnot(abs(sum(g$cells$area_m2)/as.numeric(sf::st_area(irregular))-1)<1e-7)
projected <- tempfile(fileext=".gpkg"); sf::st_write(sf::st_transform(a$aoi,3978),projected,quiet=TRUE)
cfg <- course_config;cfg$aoi_path <- projected
b <- .read_spatial_inputs(cfg)
stopifnot(nrow(b$points)==4,abs(as.numeric(sf::st_area(a$aoi))/as.numeric(sf::st_area(b$aoi))-1)<1e-6)
unlink(projected)
# Positive longitudes on the existing Canadian AOI must not silently map to Canada.
wrong <- a$cores;wrong$longitude <- abs(wrong$longitude)
f <- tempfile(fileext=".csv");write.csv(wrong,f,row.names=FALSE)
stopifnot(inherits(try(.read_spatial_inputs(course_config,f),silent=TRUE),"try-error"));unlink(f)
message("Revision checks passed: gstat equivalence, memory cap, kriging train/test provenance, failed stage cleanup, holes, projected AOI and wrong longitude signs.")
