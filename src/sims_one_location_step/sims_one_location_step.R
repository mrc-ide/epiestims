source("R/sim_utils.R")
set.seed(42)

ndays <- 100
n_loc <- 1
n_v <- 2

# SI distr
si_mu_ref <- 5.4
si_std_ref <- 1.5
si_distr_ref <- discr_si(0:30, mu = si_mu_ref, sigma = si_std_ref)
si_distr_ref <- si_distr_ref / sum(si_distr_ref)
si_no_zero_ref <- si_distr_ref[-1]
## Pierre suggests setting up the grid that
## Rt does not exceed 5 for variant or for ref.
## This is to prevent simulating unrealistically
## large numbers
sim_params <- expand.grid(
  rt_ref = c(1.4, 1.6),
  rt_post_step = c(1.1, 1.2),
  step_time = 30,
  epsilon = c(seq(from = 1, to = 2, by = 0.1), 2.5, 3),
  si_mu_variant = 1 * si_mu_ref,
  si_std_variant = si_std_ref
)

sim_params <- sim_params %>% 
  filter(!(rt_ref == 1.4 & rt_post_step == 1.2)) %>% 
  filter(!(rt_ref == 1.6 & rt_post_step == 1.1))

# Number of simulations
nsims <- ifelse(short_run, 1, 100)
rows <- ifelse(short_run, 2, nrow(sim_params))
sim_params <- sim_params[seq_len(rows), ]
incid_init <- initial_incidence()
##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
# plan(multicore)
simulated_incid <- future_pmap(
  sim_params,
  function(rt_ref, rt_post_step, step_time,
           epsilon, si_mu_variant, si_std_variant) {
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    si_no_zero_var <- si_distr_variant[-1]
    si_for_sim <- cbind(si_no_zero_ref, si_no_zero_var)
    simulate_stepwise_incid_wrapper(
      rt_ref, rt_post_step, step_time, epsilon, si_for_sim,
      incid_init = incid_init, nsims = nsims
    )
  }, .options = furrr_options(seed = TRUE)
)

si_for_est <- future_pmap(
  sim_params,
  function(rt_ref, rt_post_step, step_time,
           epsilon, si_mu_variant, si_std_variant) {
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    cbind(si_distr_ref, si_distr_variant)
  }, .options = furrr_options(seed = TRUE)
)

saveRDS(simulated_incid, "incid.rds")
saveRDS(si_for_est, "si_for_est.rds")
saveRDS(sim_params, "param_grid.rds")
