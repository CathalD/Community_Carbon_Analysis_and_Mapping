# The workshop's Gaussian conditioning/precision algebra on the original scale.
# This is an explicit response-scale adaptation for Li zero-depth/zero-carbon
# observations. No epsilon, hidden exclusion, or clipping is applied.
if (!exists("distance_matrix", mode = "function")) source("helpers/workshop_gp.R")

gaussian_results <- function(prior_mean, prior_sd, mu, variance) {
  stopifnot(all(is.finite(mu)), all(is.finite(variance)), all(variance >= -1e-8),
    all(variance <= prior_sd^2 + 1e-7*pmax(1, prior_sd^2)))
  variance <- pmax(0, variance); sd <- sqrt(variance)
  data.frame(prior_mean = prior_mean, prior_sd = prior_sd,
    posterior_mean = mu, posterior_sd = sd,
    lower95 = mu-qnorm(.975)*sd, upper95 = mu+qnorm(.975)*sd,
    variance_reduction = 1-variance/prior_sd^2,
    change_from_prior = mu-prior_mean,
    negative_mean = as.integer(mu < 0),
    probability_below_zero = ifelse(sd > 0, pnorm(-mu/sd), as.numeric(mu < 0)))
}

workshop_gaussian_gp <- function(core_xy, observed, core_mean, core_sd,
  new_xy, new_mean, new_sd, range_m = 2000, observation_sd = 0) {
  stopifnot(length(observed) > 0, all(is.finite(observed) & observed >= 0),
    all(is.finite(core_mean)), all(is.finite(core_sd) & core_sd > 0),
    all(is.finite(new_mean)), all(is.finite(new_sd) & new_sd > 0),
    range_m > 0, length(observation_sd) %in% c(1, length(observed)))
  K <- outer(core_sd, core_sd)*exp(-distance_matrix(core_xy, core_xy)/range_m)
  V <- K + diag(rep(observation_sd, length.out = length(observed))^2, nrow = length(observed))
  R <- chol(V)
  alpha <- backsolve(R, forwardsolve(t(R), observed-core_mean))
  mu <- variance <- numeric(length(new_mean))
  for (start in seq(1, length(new_mean), by = 5000)) {
    ii <- start:min(start+4999, length(new_mean))
    cross <- outer(new_sd[ii], core_sd)*exp(-distance_matrix(new_xy[ii,,drop=FALSE], core_xy)/range_m)
    mu[ii] <- new_mean[ii]+as.vector(cross %*% alpha)
    variance[ii] <- new_sd[ii]^2-colSums(forwardsolve(t(R), t(cross))^2)
  }
  gaussian_results(new_mean, new_sd, mu, variance)
}

workshop_gaussian_fusion <- function(prior_mean, prior_sd, local_mean, error_sd) {
  stopifnot(all(is.finite(local_mean)), is.finite(error_sd), error_sd > 0,
    all(is.finite(prior_sd) & prior_sd > 0))
  prior_variance <- prior_sd^2
  weight <- prior_variance/(prior_variance+error_sd^2)
  result <- gaussian_results(prior_mean, prior_sd,
    prior_mean+weight*(local_mean-prior_mean), (1-weight)*prior_variance)
  result$local_estimate <- local_mean
  result$data_weight <- weight
  result
}
