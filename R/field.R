# Read workbook inputs and calculate soil organic carbon independently of Excel.
# Density is g dry soil/cm3; carbon is percent dry mass; depths are cm.

.present <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
.number <- function(x) suppressWarnings(as.numeric(x))
.empty_qc <- function() data.frame(core_id=character(), section_id=character(),
  excel_row=integer(), check=character(), severity=character(), message=character(),
  original_value=numeric(),analysis_value=numeric())
.qc_row <- function(core_id="", section_id="", excel_row=NA_integer_, check,
                    message, severity="Review",original_value=NA_real_,analysis_value=NA_real_) {
  data.frame(core_id=as.character(core_id), section_id=as.character(section_id),
    excel_row=as.integer(excel_row), check=check, severity=severity, message=message,
    original_value=original_value,analysis_value=analysis_value)
}

validate_core_layers <- function(layers) {
  required <- c("core_id","section_id","top_cm","bottom_cm","bd_g_cm3","carbon_used_pct")
  if (!all(required %in% names(layers)) || !nrow(layers)) stop("No complete section input table supplied.")
  if (any(!.present(layers$core_id)) || any(!.present(layers$section_id))) stop("A section is missing its core or section ID.")
  if (anyDuplicated(layers$section_id)) stop("Duplicate section IDs: use a unique ID for each sample.")
  for (v in c("top_cm","bottom_cm","bd_g_cm3","carbon_used_pct")) {
    if (!is.numeric(layers[[v]]) || any(!is.finite(layers[[v]]))) stop("Missing or invalid ", v, ".")
  }
  if (any(layers$top_cm < 0 | layers$bottom_cm <= layers$top_cm)) stop("Invalid section depths.")
  if (any(layers$bd_g_cm3 <= 0)) stop("Bulk density must be positive.")
  if (any(layers$carbon_used_pct < 0 | layers$carbon_used_pct > 100)) stop("Carbon must be between 0 and 100 percent.")
  for (z in split(layers, layers$core_id)) {
    z <- z[order(z$top_cm), ]
    if (nrow(z)>1 && any(z$top_cm[-1] < head(z$bottom_cm,-1)-1e-8)) {
      stop("Overlapping sections in ", z$core_id[1], "; stocks cannot be summed safely.")
    }
  }
  invisible(TRUE)
}

complete_depth_window <- function(z, depth) {
  z <- z[order(z$top_cm), ]; z <- z[z$top_cm < depth, ]
  nrow(z)>0 && abs(z$top_cm[1])<1e-8 && max(z$bottom_cm)>=depth &&
    (nrow(z)==1 || all(abs(z$top_cm[-1]-head(z$bottom_cm,-1))<1e-8))
}

validate_example_density_override <- function(layers,override) {
  if (length(override)!=1 || (!is.na(override) && (!is.finite(override) || override<=0))) {
    stop("Invalid J7 density override.")
  }
  if (is.na(override)) return(invisible(TRUE))
  j <- which(layers$excel_row==7L)
  valid <- length(j)==1 && identical(as.character(layers$core_id[j]),"SIM-WET-01") &&
    identical(as.character(layers$section_id[j]),"SIM-WET-01-03") &&
    is.finite(layers$bd_g_cm3[j]) && any(abs(layers$bd_g_cm3[j]-c(.9,.09))<1e-10)
  if (!valid) stop("The J7 override belongs to the supplied example only (SIM-WET-01-03, original density 0.9 or 0.09). Set density_override_J7=NA_real_ in config.R for your own workbook.")
  invisible(TRUE)
}

import_core_workbook <- function(path, density_override_J7=NA_real_) {
  if (!requireNamespace("readxl", quietly=TRUE)) stop("Install readxl first.")
  if (!file.exists(path)) stop("Workbook not found: ", path)
  read_inputs <- function(sheet, columns) {
    x <- as.data.frame(readxl::read_excel(path, sheet=sheet, skip=3,
      col_types="text", .name_repair="minimal"), stringsAsFactors=FALSE)
    if (ncol(x)<columns) stop("The workbook input layout has changed in ", sheet, ".")
    x[,seq_len(columns),drop=FALSE]
  }
  soil <- read_inputs("3. Peat Data",21)
  log <- read_inputs("2. Core Log",22)
  plots <- read_inputs("1. Plot & Site Log",16)
  soil$excel_row <- seq_len(nrow(soil))+4L
  names(soil)[1:21] <- c("plot_id","core_id","drive","section_id","top_cm","bottom_cm",
    "cached_thickness","von_post","bd_basis","bd_g_cm3","water_wet_mass_pct",
    "loi_pct","oc_pct","cached_carbon_used","cached_overlap","cached_stock",
    "cached_reference_stock","cached_carbon_density","cached_cumulative","cached_qc","notes")
  names(log) <- c("core_id","plot_id","latitude","longitude","date","drives",
    "bore_depth_cm","recovered_cm","compaction_cm","recovery_ratio","organic_contact_cm",
    "contact_reached","von_post_base","chamber_length_cm","chamber_diameter_cm","photos",
    "cached_depth","cached_stock","cached_reference_stock","reference_reached","cached_qc","notes")
  names(plots) <- c("plot_id","site_id","study_area","date","location","latitude","longitude",
    "gnss_accuracy_m","datum","elevation_m","wetland_type","microform","vegetation",
    "water_table_cm","plot_area_m2","notes")
  soil <- soil[.present(soil$core_id) | .present(soil$section_id) |
    .present(soil$top_cm) | .present(soil$bottom_cm), , drop=FALSE]
  log <- log[.present(log$core_id), , drop=FALSE]
  plots <- plots[.present(plots$plot_id), , drop=FALSE]
  if (!nrow(soil) || !nrow(log) || !nrow(plots)) stop("Enter section, core and plot records first.")
  if (anyDuplicated(log$core_id) || anyDuplicated(plots$plot_id)) stop("Duplicate core or plot IDs.")
  sm <- match(soil$core_id,log$core_id); pm <- match(log$plot_id,plots$plot_id)
  if (anyNA(sm) || anyNA(pm)) stop("A section/core/plot ID has no matching record.")
  if (any(!.present(soil$plot_id) | soil$plot_id != log$plot_id[sm])) stop("Section and core plot IDs disagree.")
  for (v in c("top_cm","bottom_cm","bd_g_cm3","water_wet_mass_pct","loi_pct","oc_pct")) {
    raw <- soil[[v]]; soil[[v]] <- .number(raw)
    if (v %in% c("top_cm","bottom_cm","bd_g_cm3","oc_pct") &&
        any(.present(raw) & !is.finite(soil[[v]]))) stop("Entered ",v," must be a finite number; blank carbon may use LOI.")
  }
  soil$bd_input_g_cm3 <- soil$bd_g_cm3
  qc <- .empty_qc()
  validate_example_density_override(soil,density_override_J7)
  if (!is.na(density_override_J7)) {
    j <- which(soil$excel_row==7L)
    if (length(j)!=1) stop("The requested J7 override has no matching input row.")
    soil$bd_g_cm3[j] <- density_override_J7
    qc <- rbind(qc,.qc_row(soil$core_id[j],soil$section_id[j],7L,"density_override",
      paste0("Analysis uses confirmed J7 density ",density_override_J7,
        " g/cm3; original workbook input ",soil$bd_input_g_cm3[j]," is retained."),"Information",
      original_value=soil$bd_input_g_cm3[j],analysis_value=density_override_J7))
  }
  settings <- readxl::read_excel(path,sheet="7. Settings",range="A5:B18",
    col_names=FALSE,col_types="text",.name_repair="minimal")
  factor <- .number(settings[[2]][match("CARBON_FRACTION_OM",settings[[1]])])
  missing_oc <- !is.finite(soil$oc_pct)
  if (any(missing_oc) && (length(factor)!=1 || !is.finite(factor) || factor<=0 || factor>1)) {
    stop("A valid LOI-to-carbon factor is required when organic carbon is absent.")
  }
  soil$carbon_basis <- ifelse(missing_oc,"LOI-derived organic carbon","Entered organic carbon")
  soil$carbon_used_pct <- soil$oc_pct
  soil$carbon_used_pct[missing_oc] <- soil$loi_pct[missing_oc]*factor
  validate_core_layers(soil)
  for (v in c("latitude","longitude")) {
    raw <- log[[v]]; log[[v]] <- .number(raw)
    if (any(.present(raw) & !is.finite(log[[v]]))) stop("Entered ",v," must be a finite number.")
  }
  plots$latitude <- .number(plots$latitude); plots$longitude <- .number(plots$longitude)
  if (any(is.finite(log$latitude) & abs(log$latitude)>90) ||
      any(is.finite(log$longitude) & abs(log$longitude)>180)) stop("Coordinates are outside latitude/longitude bounds.")
  log$site_id <- plots$site_id[pm]
  log$wetland_type <- plots$wetland_type[pm]
  log$vegetation <- plots$vegetation[pm]
  log$study_area <- plots$study_area[pm]
  # Explicit site labels supplied with the example, independent of core IDs.
  site_ecosystems <- c("SIM-BOG"="Bog","SIM-FEN"="Fen",
    "SIM-FOREST-A"="Upland conifer forest","SIM-FOREST-B"="Upland mixed forest")
  log$ecosystem <- unname(site_ecosystems[log$site_id])
  habitat <- trimws(tolower(ifelse(is.na(log$wetland_type),"",log$wetland_type)))
  log$ecosystem[habitat=="bog"] <- "Bog"
  log$ecosystem[habitat=="fen"] <- "Fen"
  explicit_forest <- habitat %in% c("forest","upland forest","conifer forest","mixed forest")
  log$ecosystem[explicit_forest] <- log$wetland_type[explicit_forest]
  log$ecosystem[is.na(log$ecosystem)] <- "Not recorded"
  log$group <- ifelse(log$ecosystem %in% c("Bog","Fen"),"Wetland cores",
    ifelse(grepl("forest",log$ecosystem,ignore.case=TRUE),"Forest cores","Other cores"))
  soil$ecosystem <- log$ecosystem[match(soil$core_id,log$core_id)]
  soil$group <- log$group[match(soil$core_id,log$core_id)]
  for (i in seq_len(nrow(soil))) {
    z <- soil[i,]
    add <- function(check,message,severity="Review") .qc_row(z$core_id,z$section_id,z$excel_row,check,message,severity)
    if (is.finite(z$loi_pct) && z$carbon_used_pct > z$loi_pct+.01) qc <- rbind(qc,add("carbon_exceeds_loi",
      "Entered organic carbon exceeds LOI. Organic carbon is used; review the independent LOI value."))
    if (is.finite(z$water_wet_mass_pct)) {
      w <- z$water_wet_mass_pct
      if (w<0 || w>=100 || z$bd_g_cm3*w/(100-w)>1+1e-8) qc <- rbind(qc,add("water_density_inconsistent",
        "Wet-mass water percentage and dry bulk density imply more water than sample volume, or water is out of range. Stocks use density and carbon only."))
    }
    if (missing_oc[i]) qc <- rbind(qc,add("loi_derived_carbon",
      paste0("Organic carbon was absent; LOI multiplied by ",factor," was used."),"Information"))
    # Notes are retained verbatim; do not silently update old nitrogen or C:N.
    note <- ifelse(is.na(z$notes),"",z$notes)
    parse_note <- function(pattern) {
      m <- regexec(pattern,note); x <- regmatches(note,m)[[1]]
      if (length(x)>1) .number(x[2]) else NA_real_
    }
    n_pct <- parse_note("N%=([0-9.]+)"); cn <- parse_note("C:N=([0-9.]+)")
    if (is.finite(n_pct) && n_pct>0 && is.finite(cn) &&
        abs(z$carbon_used_pct/n_pct-cn)>max(.2,abs(cn)*.02)) qc <- rbind(qc,add("notes_cn_inconsistent",
          "C:N in the retained note does not agree with entered carbon and the note's nitrogen percentage."))
  }
  for (z in split(soil,soil$core_id)) {
    if (!complete_depth_window(z,max(z$bottom_cm))) qc <- rbind(qc,.qc_row(z$core_id[1],
      check="profile_gaps",message="Profile has missing surface coverage or gaps. Sampled stocks sum entered sections; fixed-depth stocks require complete coverage."))
  }
  for (i in which(!is.finite(log$latitude) | !is.finite(log$longitude))) qc <- rbind(qc,.qc_row(log$core_id[i],
    check="missing_coordinates",message="Coordinates are missing; nonspatial stock calculation remains available."))
  labels <- c(plots$study_area,plots$location,plots$notes,log$notes,soil$notes)
  if (any(grepl("simulat|synthetic|hypothetical|unverified",labels,ignore.case=TRUE),na.rm=TRUE)) {
    qc <- rbind(qc,.qc_row(check="provenance_labels",
      message="Workbook retains simulated or hypothetical labels and locations. Use the configured data origin; the labels do not establish measured field provenance.",severity="Information"))
  }
  soil <- soil[order(soil$core_id,soil$top_cm), , drop=FALSE]
  rownames(soil) <- NULL; rownames(log) <- NULL
  list(layers=soil,cores=log,plots=plots,qc=qc,source=normalizePath(path))
}

calculate_core_stocks <- function(layers,cores) {
  validate_core_layers(layers)
  if (anyDuplicated(cores$core_id)) stop("Duplicate core metadata IDs.")
  if (anyNA(match(layers$core_id,cores$core_id))) stop("Missing core metadata.")
  layers <- layers[order(layers$core_id,layers$top_cm), , drop=FALSE]
  layers$thickness_cm <- layers$bottom_cm-layers$top_cm
  # 1 g/cm2 = 10 kg/m2; percent is divided by 100 exactly once.
  layers$stock_kg_m2 <- with(layers,bd_g_cm3*carbon_used_pct/100*thickness_cm*10)
  layers$cumulative_sampled_stock_kg_m2 <- ave(layers$stock_kg_m2,layers$core_id,FUN=cumsum)
  result <- do.call(rbind,lapply(split(layers,layers$core_id),function(z) {
    m <- cores[match(z$core_id[1],cores$core_id), , drop=FALSE]
    window_stock <- function(depth) {
      if (!complete_depth_window(z,depth)) return(NA_real_)
      sum(z$bd_g_cm3*z$carbon_used_pct/100*pmax(0,pmin(z$bottom_cm,depth)-z$top_cm)*10)
    }
    data.frame(core_id=m$core_id,plot_id=m$plot_id,ecosystem=m$ecosystem,group=m$group,
      latitude=m$latitude,longitude=m$longitude,sampled_depth_cm=max(z$bottom_cm),n_sections=nrow(z),
      complete_sampled_profile=complete_depth_window(z,max(z$bottom_cm)),
      stock_0_30_kg_m2=window_stock(30),stock_0_100_kg_m2=window_stock(100),
      sampled_stock_kg_m2=sum(z$stock_kg_m2),complete_to_100cm=complete_depth_window(z,100))
  }))
  result$stock_0_30_t_ha <- result$stock_0_30_kg_m2*10
  rownames(result) <- NULL
  list(cores=result,layers=layers)
}

analyse_core_data <- function(input,out="outputs/part1",data_origin="User-supplied") {
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 first.")
  results <- calculate_core_stocks(input$layers,input$cores)
  layers <- results$layers; cores <- results$cores
  layers$data_origin <- data_origin; cores$data_origin <- data_origin
  statistics <- do.call(rbind,lapply(c("group","ecosystem"),function(level) {
    do.call(rbind,lapply(split(cores,cores[[level]]),function(z) {
      x <- z$stock_0_30_kg_m2; x <- x[is.finite(x)]
      avg <- if(length(x)) mean(x) else NA_real_; s <- if(length(x)>1) stats::sd(x) else NA_real_
      data.frame(aggregation=level,group=z[[level]][1],carbon_pool="Soil organic carbon",depth_cm="0-30",
        n_cores=length(x),mean_kg_m2=avg,sample_sd_kg_m2=s,
        cv_pct=if(is.finite(avg) && avg!=0) s/avg*100 else NA_real_,
        mean_t_ha=avg*10,sample_sd_t_ha=s*10)
    }))
  }))
  statistics$data_origin <- data_origin
  qc <- input$qc
  for (i in seq_len(nrow(cores))) {
    old <- .number(input$cores$cached_stock[match(cores$core_id[i],input$cores$core_id)])
    if (length(old)==1 && is.finite(old) && abs(cores$sampled_stock_kg_m2[i]-old)>1e-6) {
      qc <- rbind(qc,.qc_row(cores$core_id[i],check="cached_stock_differs",
        message="Independent stock differs from the workbook cached total. Review configured overrides or recalculate Excel.",severity="Information"))
    }
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for (item in c("layer_stocks","core_stocks","descriptive_statistics","qc_flags")) {
    value <- switch(item,layer_stocks=layers,core_stocks=cores,descriptive_statistics=statistics,qc_flags=qc)
    utils::write.csv(value,file.path(out,paste0(item,".csv")),row.names=FALSE,na="")
  }
  theme <- ggplot2::theme_bw(base_size=12)
  profile <- rbind(data.frame(core_id=layers$core_id,top=layers$top_cm,bottom=layers$bottom_cm,
      variable="Organic C (%)",value=layers$carbon_used_pct),
    data.frame(core_id=layers$core_id,top=layers$top_cm,bottom=layers$bottom_cm,
      variable="Bulk density (g/cm³)",value=layers$bd_g_cm3))
  profile$depth_cm <- (profile$top+profile$bottom)/2
  # Break profile lines at gaps. Endpoints extend only within sampled sections.
  lines <- do.call(rbind,lapply(split(profile,interaction(profile$core_id,profile$variable,drop=TRUE)),function(z) {
    z <- z[order(z$top), ]; z$segment <- cumsum(c(TRUE,abs(z$top[-1]-head(z$bottom,-1))>1e-8))
    do.call(rbind,lapply(split(z,z$segment),function(y) {
      start <- y[1, ]; start$depth_cm <- start$top
      end <- y[nrow(y), ]; end$depth_cm <- end$bottom
      rbind(start,y,end)
    }))
  }))
  p <- ggplot2::ggplot(profile,ggplot2::aes(value,depth_cm,colour=core_id))+
    ggplot2::geom_path(data=lines,ggplot2::aes(group=interaction(core_id,segment)),linewidth=.8)+
    ggplot2::geom_point(size=2)+ggplot2::scale_y_reverse(limits=c(max(layers$bottom_cm),0))+
    ggplot2::facet_wrap(~variable,scales="free_x")+theme+
    ggplot2::labs(x=NULL,y="Depth (cm)",colour="Core")
  ggplot2::ggsave(file.path(out,"depth_profiles.png"),p,width=10,height=6,dpi=160)
  curve <- do.call(rbind,lapply(split(layers,layers$core_id),function(z) {
    z <- z[order(z$top_cm), ]; z$segment <- cumsum(c(TRUE,abs(z$top_cm[-1]-head(z$bottom_cm,-1))>1e-8))
    do.call(rbind,lapply(split(z,z$segment),function(y) {
      start <- y[1, ]; start$bottom_cm <- start$top_cm
      start$cumulative_sampled_stock_kg_m2 <- start$cumulative_sampled_stock_kg_m2-start$stock_kg_m2
      rbind(start,y)
    }))
  }))
  p <- ggplot2::ggplot(curve,ggplot2::aes(cumulative_sampled_stock_kg_m2,bottom_cm,colour=core_id,
      group=interaction(core_id,segment)))+ggplot2::geom_path(linewidth=.9)+ggplot2::geom_point(size=1.8)+
    ggplot2::scale_y_reverse()+theme+ggplot2::labs(x="kg C/m²",
      y="Depth (cm)",colour="Core")
  ggplot2::ggsave(file.path(out,"cumulative_stock.png"),p,width=9,height=6,dpi=160)
  p <- ggplot2::ggplot(cores,ggplot2::aes(core_id,stock_0_30_kg_m2,fill=group))+
    ggplot2::geom_col(width=.65,na.rm=TRUE)+theme+
    ggplot2::labs(x=NULL,y="0–30 cm (kg C/m²)",fill=NULL)
  ggplot2::ggsave(file.path(out,"common_depth_stock.png"),p,width=9,height=5,dpi=160)
  list(cores=cores,layers=layers,statistics=statistics,qc=qc,metadata=input$cores)
}
