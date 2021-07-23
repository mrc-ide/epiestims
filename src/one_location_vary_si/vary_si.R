## orderly::orderly_develop_start(use_draft = "newer", parameters = list(short_run = TRUE))
simulated_incid <- readRDS("incid.rds")
si_for_est <- readRDS("si_for_est.rds")

priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 15000L, burnin = 7500L, thin = 20L
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
      attempt <- 1
      while (! out[["convergence"]]) {
        message("Attempt ", attempt)
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
        attempt <- attempt + 1
        ## so that we don't end in an infinite loop
        if (attempt > 3) {
          message("Aborting after 3 attempts")
          break
        }
      }
      list(out, out[["convergence"]])
    }
    )
    }
    )
  }
)

saveRDS(results, "estimate_joint_output.rds")

