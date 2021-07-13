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

## Seed with 1 case but select simulations
## that went on to generate at least 20 cases
## Ref and variant have different initial incidence conditions, saved in this list
initial_incidence <- list(
  incidence::incidence(rep(seq(1, 10), each = 20)),
  incidence::incidence(rep(1, 1))
)

sim_params <- expand.grid(
  rt_ref = c(1.2, 3),
  epsilon = c(seq(from = 1, to = 2, by = 0.1), 2.5, 3),
  si_mu_variant = c(0.5, 0.75, 1, 1.25, 1.5) * si_mu_ref,
  si_std_variant = si_std_ref
)

##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
simulated_incid <- pmap(
  sim_params,
  function(rt_ref, epsilon, si_mu_variant, si_std_variant) {
    ## Calculate reproduction number for variant
    rt_variant <- epsilon * rt_ref
    ## Assume reproduction number remains the same
    ## over the time period
    ## Make a vector that goes across rows
    R <- array(NA, dim = c(ndays, n_loc, n_v))
    R[,,1] <- rep(rt_ref, each = ndays)
    R[,,2] <- rep(rt_variant, each = ndays)
    si_distr_variant <- discr_si(
      0:30, mu = si_mu_variant, sigma = si_std_variant
    )
    si_distr_variant <- si_distr_variant / sum(si_distr_variant)
    si_no_zero_var <- si_distr_variant[-1]
    si_for_sim <- cbind(si_no_zero_ref, si_no_zero_var)
    ## Because we are starting with a small seed
    ## we simulate 10 times as many trajectories
    ## as we need so that we have nsim after
    ## filtering
    out <- imap(
      seq_len(10 * nsims), function(x, index) {
        message("Simulation ", index)
        simulate_incidence(
          initial_incidence, n_loc, n_v, ndays, R, si_for_sim
        )
      }
    )
    ## total number of cases at the end of the
    ## simulation for the variant.
    ncases1 <- map_dbl(out, function(x) sum(x[, , 1]))
    ncases2 <- map_dbl(out, function(x) sum(x[, , 2]))
    out <- out[ncases1 > 20 & ncases2 > 20]
    message("# of simulations with more than 20 cases ", length(out))
    ## At this point out will have either less than
    ## or more than the desired number of simulations
    ## Sample to make sure we have exactly nsim
    index <- sample(length(out), size = nsims, replace = TRUE)
    out <- out[index]
    names(out) <- seq_len(nsims)
    out
  }
)


tmax_all <- seq(20, 60, by = 10)
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
      EpiEstim:::estimate_joint(
        x, si_for_est, priors, seed = 1,
        t_min = 15L, t_max = as.integer(tmax),
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
