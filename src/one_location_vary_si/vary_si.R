source("R/sim_utils.R")
set.seed(42)

ndays <- 100
n_loc <- 1
n_v <- 2

# SI distr

si_mu_ref <- 6.83
si_std_ref <- 3.8
si_distr_ref <- discr_si(0:30, mu = si_mu_ref, sigma = si_std_ref)
si_distr_ref <- si_distr_ref / sum(si_distr_ref)
si_no_zero_ref <- si_distr_ref[-1]

# Number of simulations
nsims <- ifelse(short_run, 1, 100)

## Other common things
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 20000L, burnin = 1000L, thin = 40L
)


sim_params <- expand.grid(
  rt_ref = 3,
  ##epsilon = c(seq(from = 1, to = 2, by = 0.1), 2.5, 3),
  epsilon = 1,
  si_mu_variant = 1.5 * si_mu_ref,
  si_std_variant = si_std_ref
)

rows <- ifelse(short_run, 2, nrow(sim_params))
sim_params <- sim_params[seq_len(rows), ]
incid_init <- initial_incidence()
##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
simulated_incid <- pmap(
  sim_params,
  function(rt_ref, epsilon, si_mu_variant, si_std_variant) {
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    si_no_zero_var <- si_distr_variant[-1]
    si_for_sim <- cbind(si_no_zero_ref, si_no_zero_var)
    simulate_incid_wrapper(
      rt_ref, epsilon, si_for_sim, incid_init = incid_init, nsims = nsims)
  }
)


tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all

## Estimate epsilon
results <- pmap(
  list(
    incid = simulated_incid,
    si_mu_variant = sim_params$si_mu_variant,
    si_std_variant = sim_params$si_std_variant
  ),
  function(incid, si_mu_variant, si_std_variant) {
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    si_for_est <- cbind(si_distr_ref, si_distr_variant)
    map(tmax_all, function(tmax) {
    message("tmax = ", tmax)
    ## Loop over the first dimension which is
    ## the set of simulations
    map(incid, function(x) {
      t_min <- EpiEstim::compute_t_min(x, si_for_est)
      t_max <- as.integer(t_min + tmax)
      t_max <- min(t_max, nrow(x))
      EpiEstim:::estimate_joint(
        x, si_for_est, priors, seed = 1,
        t_min = t_min,
        t_max = t_max,
        mcmc_control = mcmc_controls
      )

    }
    )
    }
    )
  }
)


saveRDS(sim_params, "param_grid.rds")
saveRDS(results, "estimate_joint_output.rds")
saveRDS(simulated_incid, "sim_incid.rds")
