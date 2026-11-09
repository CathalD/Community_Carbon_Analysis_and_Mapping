# Fixed exploratory settings. No selection using community outcomes.
# Environment: exp(-0.5 * sum((z_target - z_donor)^2 / ell^2)).
# z uses donor-only median/MAD; ell=sqrt(5) for the five prespecified predictors.
environment_lengthscale <- sqrt(5)
top_donors_to_export <- 25L
# AlphaEarth bandwidth is the median nearest-neighbour chord distance in each
# training pool. All 64 vector dimensions are L2-normalized without test fitting.
# No measurement-error or predictive variance is assigned to donor spread.
