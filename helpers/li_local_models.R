# Local models from workshop Scripts 2/3, applied to Li-matched organic variables.
# This file uses measured layers and explicit contacts, never Route L model output.

li_prepare_field_observations <- function(layers, locations, contacts,
    mode = "all_eight_layer_evidence", organic_soc_threshold_g_kg = 80) {
  stopifnot(mode %in% c("all_eight_layer_evidence", "confirmed_only"),
    !anyDuplicated(locations$plot_id), !anyDuplicated(contacts$plot_id),
    all(is.finite(layers$soc) & layers$soc >= 0 & layers$soc <= 1000),
    all(is.finite(layers$bulk_density) & layers$bulk_density > 0),
    all(is.finite(layers$coarse_frag) & layers$coarse_frag >= 0 & layers$coarse_frag <= 1),
    all(layers$depth_from >= 0 & layers$depth_to > layers$depth_from))
  out <- lapply(seq_len(nrow(locations)), function(i) {
    id <- as.character(locations$plot_id[i])
    d <- layers[layers$plot_id == id, , drop = FALSE]
    if (!nrow(d)) stop("No measured layers for ", id)
    d <- d[order(d$depth_from, d$depth_to), , drop = FALSE]
    if (abs(d$depth_from[1]) > 1e-8 || (nrow(d) > 1 &&
        any(abs(d$depth_from[-1] - head(d$depth_to, -1)) > 1e-8))) {
      stop("Organic integration needs continuous surface coverage: ", id)
    }
    j <- match(id, contacts$plot_id)
    confirmed <- !is.na(j)
    has_mineral <- any(d$soc < organic_soc_threshold_g_kg)
    if (confirmed) {
      depth <- contacts$organic_horizon_depth_cm[j]
      if (!is.finite(depth) || depth < 0 || depth > max(d$depth_to) + 1e-8) {
        stop("Confirmed contact is invalid or beyond measured support at ", id)
      }
      evidence <- "field_confirmed_contact"
      detail <- as.character(contacts$contact_method[j])
      complete_organic_profile <- TRUE
    } else if (has_mineral) {
      depth <- d$depth_from[which(d$soc < organic_soc_threshold_g_kg)[1]]
      evidence <- "layer_classified_contact"
      detail <- paste0("First layer with SOC < ", organic_soc_threshold_g_kg,
        " g/kg; boundary at its top; classification-derived, not field-confirmed")
      complete_organic_profile <- TRUE
    } else {
      depth <- max(d$depth_to)
      evidence <- "organic_to_sampled_bottom_lower_bound"
      detail <- "No mineral layer or confirmed contact; sampled bottom is a lower bound"
      complete_organic_profile <- FALSE
    }
    stock <- sum(d$soc * d$bulk_density * (1 - d$coarse_frag) *
      pmax(0, pmin(d$depth_to, depth) - d$depth_from) / 100)
    canfit <- complete_organic_profile && (confirmed || mode == "all_eight_layer_evidence")
    data.frame(core_id = id, longitude = locations$longitude[i], latitude = locations$latitude[i],
      depth = depth, carbon = stock, peat_depth_cm = depth, organic_carbon_kg_m2 = stock,
      evidence = evidence, confirmed = confirmed, canfit = canfit,
      complete_organic_profile = complete_organic_profile,
      zero_organic_depth = depth == 0, sampled_bottom_cm = max(d$depth_to),
      contact_rule = detail, contact_mode = mode, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

li_cv_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  if (!any(keep)) return(c(RMSE = NA_real_, MAE = NA_real_, bias = NA_real_, R2 = NA_real_))
  observed <- observed[keep]; predicted <- predicted[keep]
  errors <- predicted - observed
  ss <- sum((observed - mean(observed))^2)
  c(RMSE = sqrt(mean(errors^2)), MAE = mean(abs(errors)), bias = mean(errors),
    R2 = if (ss > 0) 1-sum(errors^2)/ss else NA_real_)
}

li_idw_predict <- function(training_points, prediction_points, idp, nmax) {
  as.numeric(gstat::idw(response ~ 1, locations = training_points,
    newdata = prediction_points, idp = idp, nmax = min(nmax, nrow(training_points)),
    debug.level = 0)$var1.pred)
}

li_idw_models <- function(field, core_points, grid_points) {
  n <- nrow(field)
  if (n < 2) stop("IDW needs at least two local observations for a descriptive holdout")
  candidates <- expand.grid(idp = c(.5, 1, 1.5, 2, 3, 4), nmax = 2:5)
  candidates <- candidates[candidates$nmax <= n-1L, , drop = FALSE]
  cv <- function(idp, nmax) vapply(seq_len(n), function(i) {
    li_idw_predict(core_points[-i, ], core_points[i, ], idp, nmax)
  }, numeric(1))
  fallback <- nrow(candidates) == 0
  if (fallback) {
    candidates <- data.frame(idp = 1, nmax = min(2, n), RMSE = NA_real_)
    best <- candidates
    heldout <- cv(best$idp, best$nmax)
    best$RMSE <- sqrt(mean((heldout-field$response)^2))
    status <- "Fixed IDW fallback: standard workshop tuning unavailable with fewer than three cores"
  } else {
    candidates$RMSE <- vapply(seq_len(nrow(candidates)), function(i) {
      sqrt(mean((cv(candidates$idp[i], candidates$nmax[i])-field$response)^2))
    }, numeric(1))
    best <- candidates[which.min(candidates$RMSE), , drop = FALSE]
    heldout <- cv(best$idp, best$nmax)
    status <- "Workshop LOOCV tuning; tuning and reported error use the same cores (not nested)"
  }
  best$N_Cores <- n
  best$fallback <- fallback
  best$status <- status
  list(parameters = best, tuning = candidates, heldout = heldout,
    mapped = li_idw_predict(core_points, grid_points, best$idp, best$nmax), status = status)
}

li_regression_models <- function(model_data, candidates) {
  availability <- list(); fits <- list(); heldouts <- list(); evaluations <- list()
  n <- nrow(model_data)
  for (name in names(candidates)) {
    vars <- candidates[[name]]
    available <- all(vars %in% names(model_data))
    reason <- if (available) "Eligible" else "Predictor absent from prepared data"
    if (available && length(vars)) {
      available <- all(vapply(vars, function(v) all(is.finite(model_data[[v]])) &&
        length(unique(model_data[[v]])) >= 3, logical(1)))
      if (!available) reason <- "Predictor incomplete or fewer than three distinct values; all cores retained"
    }
    if (available && length(vars) > 0 && n < length(vars)+3) {
      available <- FALSE
      reason <- "Too few cores for workshop regression and whole-core holdouts"
    }
    formula <- reformulate(if (length(vars)) vars else "1", response = "response")
    if (available) {
      fit <- lm(formula, data = model_data, na.action = na.fail)
      predicted <- rep(NA_real_, n)
      for (i in seq_len(n)) {
        folded <- lm(formula, data = model_data[-i, , drop = FALSE], na.action = na.fail)
        if (anyNA(coef(folded))) {
          available <- FALSE
          reason <- "A held-out fit is rank-deficient"
          break
        }
        predicted[i] <- as.numeric(predict(folded, model_data[i, , drop = FALSE]))
      }
      if (available) {
        fits[[name]] <- fit
        heldouts[[name]] <- predicted
        K <- attr(logLik(fit), "df")
        aic <- AIC(fit)
        aicc <- if (n > K+1 && is.finite(aic)) aic+2*K*(K+1)/(n-K-1) else NA_real_
        score <- li_cv_metrics(model_data$response, predicted)
        evaluations[[name]] <- data.frame(model = name, n_cores = n,
          RMSE = score["RMSE"], MAE = score["MAE"], bias = score["bias"], R2 = score["R2"],
          parameters = K, AIC = aic, AICc = aicc, row.names = NULL,
          status = if (n < 4) "Descriptive small-n fallback; AICc unavailable" else
            "Workshop fixed-model holdouts; selection not nested")
      }
    }
    availability[[name]] <- data.frame(model = name, predictors = paste(vars, collapse = " + "),
      included = available, reason = reason, stringsAsFactors = FALSE)
  }
  if (!length(fits)) stop("No estimable regression candidate, including mean-only baseline")
  evaluation <- do.call(rbind, evaluations); rownames(evaluation) <- NULL
  finite <- which(is.finite(evaluation$AICc))
  fallback <- length(finite) == 0
  if (fallback) {
    if (!"Mean_only" %in% names(fits)) stop("AICc unavailable and mean-only fallback absent")
    selected <- "Mean_only"
    rule <- "Mean-only fallback; insufficient n for AICc selection"
  } else {
    ranking <- finite[order(evaluation$AICc[finite], evaluation$RMSE[finite])]
    selected <- evaluation$model[ranking[1]]
    rule <- "Lowest AICc, as in workshop Script 3"
  }
  evaluation$delta_AICc <- if (length(finite)) evaluation$AICc-min(evaluation$AICc[finite]) else NA_real_
  evaluation <- evaluation[order(evaluation$AICc, evaluation$RMSE, na.last = TRUE), , drop = FALSE]
  selection <- evaluation[evaluation$model == selected, , drop = FALSE]
  selection$selection_rule <- rule
  selection$fallback <- fallback
  list(fits = fits, heldouts = heldouts, evaluation = evaluation,
    availability = do.call(rbind, availability), selected = selected,
    selection = selection, fit = fits[[selected]], heldout = heldouts[[selected]])
}
