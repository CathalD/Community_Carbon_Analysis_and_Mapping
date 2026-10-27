# Optional maintainer refresh; not called by the course.
library(sf)
files <- c(countries="ne_50m_admin_0_countries",lakes="ne_50m_lakes",
  marine="ne_50m_geography_marine_polys",places="ne_50m_populated_places_simple")
cache <- tempfile("natural-earth-"); dir.create(cache)
checks <- list()
old <- sf_use_s2(FALSE)
for(layer in names(files)) {
  name <- paste0(files[[layer]],".geojson")
  path <- file.path(cache,name)
  url <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/",name)
  download.file(url,path,mode="wb")
  z <- st_read(path,quiet=TRUE)
  z <- suppressWarnings(st_crop(st_make_valid(z),c(xmin=-170,ymin=35,xmax=-40,ymax=88)))
  z <- z[,if(layer=="countries") "ADMIN" else "name",drop=FALSE]
  st_write(z,"data/context/north_america.gpkg",layer=layer,delete_layer=TRUE,quiet=TRUE)
  checks[[layer]] <- data.frame(source=name,url=url,md5=unname(tools::md5sum(path)))
}
sf_use_s2(old)
write.csv(do.call(rbind,checks),"data/context/source_checksums.csv",row.names=FALSE)
unlink(cache,recursive=TRUE)
