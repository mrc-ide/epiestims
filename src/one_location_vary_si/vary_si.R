simulated_incid <- readRDS("incid.rds")
si_for_est <- readRDS("si_for_est.rds")

priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 5000L, burnin = 2500L, thin = 10L
)

tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all

## Estimate epsilon
results <- map2(
  simulated_incid, si_for_est,
  function(incid, si) {
    map(tmax_all, function(tmax) {
    message("tmax = ", tmax)
    ## Loop over the first dimension which is
    ## the set of simulations
    map(incid, function(x) {
      t_min <- EpiEstim::compute_t_min(x, si)
      t_max <- as.integer(t_min + tmax)
      t_max <- min(t_max, nrow(x))
      out <- EpiEstim:::estimate_joint(
        x, si, priors, seed = 1,
        t_min = t_min,
        t_max = t_max,
        mcmc_control = mcmc_controls
        )
      attempts <- 1
      while (! out$convergence) {
        message("Not yet converged")
        mcmc_controls <- lapply(
          mcmc_controls, function(x) x * 2L
        )
        out <- EpiEstim:::estimate_joint(
          x, si, priors, seed = 1,
          t_min = t_min,
          t_max = t_max,
          mcmc_control = mcmc_controls
          )
        attempts <- attempts + 1
        ## so that we don't end in an infinite loop
        if (attempts > 3) {
          message("Aborting after 3 attempts")
          break
        }
      }
      out
    }
    )
    }
    )
  }
)

saveRDS(results, "estimate_joint_output.rds")

