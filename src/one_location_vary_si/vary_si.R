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
  n_iter = 5000L,
  burnin = as.integer(floor(5e3 / 2)), # speed up
  thin = 10L
)


sim_params <- expand.grid(
  rt_ref = c(1.2, 3),
  epsilon = c(seq(from = 1, to = 2, by = 0.1), 2.5, 3),
  si_mu_variant = c(0.5, 0.75, 1, 1.25, 1.5) * si_mu_ref,
  si_std_variant = si_std_ref
)
sim_params <- ifelse(
  short_run, sim_params[1:2, ], sim_params
)
incid_init <- initial_incidence()
##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
simulated_incid <- clusterMap(
  NULL,
  function(rt_ref, epsilon, si_mu_variant, si_std_variant) {
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    si_no_zero_var <- si_distr_variant[-1]
    si_for_sim <- cbind(si_no_zero_ref, si_no_zero_var)
    simulate_incid_wrapper(
      rt_ref, epsilon, si_for_sim, incid_init = incid_init, nsims = nsims)
  },
  sim_params$rt_ref, sim_params$epsilon,
  sim_params$si_mu_variant,
  sim_params$si_std_variant
)


tmax_all <- seq(20, 60, by = 10)
names(tmax_all) <- tmax_all

## Estimate epsilon
results <- clusterMap(
  NULL,
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
      EpiEstim:::estimate_joint(
        x, si_for_est, priors, seed = 1,
        t_min = NULL
        t_max = tmax,
        mcmc_control = mcmc_controls
      )
    }
    )
    }
    )
  },
  simulated_incid, sim_params$si_mu_variant, sim_params$si_std_variant
)


saveRDS(sim_params, "param_grid.rds")
saveRDS(results, "estimate_joint_output.rds")
saveRDS(simulated_incid, "sim_incid.rds")
