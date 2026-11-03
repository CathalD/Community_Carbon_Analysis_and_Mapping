# The exponential log-Gaussian update from workshop Script 5, factored for Li variables.
# mu=log(mean)-log1p((sd/mean)^2)/2; Kij=si*sj*exp(-distance/range).
# No Route L depth blend or carbon curve/tail implementation is used.
distance_matrix <- function(a, b) sqrt(outer(a[,1], b[,1], "-")^2 + outer(a[,2], b[,2], "-")^2)
log_parameters <- function(mean, sd) {
  stopifnot(all(is.finite(mean) & mean > 0), all(is.finite(sd) & sd > 0))
  variance <- log1p((sd / mean)^2)
  list(mu = log(mean) - variance/2, variance = variance)
}
workshop_gp <- function(core_xy, observed, core_mean, core_sd, new_xy, new_mean, new_sd,
                        range_m = 2000, observation_log_sd = 0) {
  stopifnot(length(observed) > 0, all(is.finite(observed) & observed > 0),
    range_m > 0, length(observation_log_sd) %in% c(1, length(observed)))
  pc <- log_parameters(core_mean, core_sd)
  pn <- log_parameters(new_mean, new_sd)
  K <- outer(sqrt(pc$variance), sqrt(pc$variance)) * exp(-distance_matrix(core_xy, core_xy)/range_m)
  V <- K + diag(rep(observation_log_sd, length.out = length(observed))^2, nrow = length(observed))
  R <- chol(V)
  alpha <- backsolve(R, forwardsolve(t(R), log(observed) - pc$mu))
  mu <- variance <- numeric(length(new_mean))
  for (start in seq(1, length(new_mean), by = 5000)) {
    ii <- start:min(start+4999, length(new_mean))
    cross <- outer(sqrt(pn$variance[ii]), sqrt(pc$variance)) *
      exp(-distance_matrix(new_xy[ii,,drop=FALSE], core_xy)/range_m)
    mu[ii] <- pn$mu[ii] + as.vector(cross %*% alpha)
    v <- pn$variance[ii] - colSums(forwardsolve(t(R), t(cross))^2)
    if (any(v < -1e-8*pmax(1, pn$variance[ii]))) stop("Invalid posterior variance")
    variance[ii] <- pmax(0, v)
  }
  data.frame(prior_mean = new_mean, prior_sd = new_sd,
    posterior_mean = exp(mu+variance/2), posterior_median = exp(mu),
    posterior_sd = sqrt(expm1(variance)*exp(2*mu+variance)),
    lower95 = exp(mu-qnorm(.975)*sqrt(variance)), upper95 = exp(mu+qnorm(.975)*sqrt(variance)),
    log_variance_reduction = 1-variance/pn$variance,
    change_from_prior = exp(mu+variance/2)-new_mean)
}
