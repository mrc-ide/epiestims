# Scenario 1: 
# 1 loc, 2 var, constant R
require(EpiEstim)
require(incidence)
require(projections)
require(purrr)
require(ggplot2)
source("simulation_functions.R")

n_v <- 2 # 2 variants
n_loc <- 1 # 1 location
ndays <- 100 # 100 time steps

## Number of simulations for each parameter
nsims <- 10

## Define range of tmax values (10 to 60) to explore
tmax_all <- seq(10, 60, 10)
tmax_all <- as.integer(tmax_all)
names(tmax_all) <- tmax_all

# Rt ref either 1.2 or 3
rt_ref <- c(1.2)

# Range of transmission advantage values to explore
transmission_advantage <- c(seq(1,2,0.1),2.5,3.0)
names(transmission_advantage) <- transmission_advantage

# SI distribution (Covid SI from IBM)
si_mean <- 6.83
si_std <- 3.8
si <- discr_si(0:30, mu = si_mean, sigma = si_std)
si <- si / sum(si)
si_no_zero <- si[-1]

# SI distribution for 2 variants
si_distr <- cbind(si_no_zero, si_no_zero)

# For estimation
si_est <- cbind(si,si)
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 10000L, burnin = as.integer(floor(1e4 / 2)),
  thin = 10L
)


# Simulate incidence 100 times

simulated_incid <- map(
  transmission_advantage, function(advantage) {
    message("transmission advantage = ", advantage)
    ## Calculate reproduction number for variant
    rt_variant <- advantage * rt_ref
    ## Constant reproduction number over the time period
    R <- array(NA, dim = c(ndays, n_loc, n_v))
    R[,,1] <- rep(rt_ref, each = ndays)
    R[,,2] <- rep(rt_variant, each = ndays)
    ##############################################################################
    ## Simulate epidemic incidence data with input reproduction numbers and si  ##
    ##############################################################################
    ## start with 20 infected individuals
    initial_incidence <- incidence::incidence(rep(1, 20))
    simulate_incidence(
      initial_incidence, n_loc, n_v, 
      ndays, R, si_distr, nsims=nsims
    )
  })


str(simulated_incid)


## Results

results <- imap(
  simulated_incid, function(incid) {
    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)
      ## Loop over the first dimension which is
      ## the set of simulations
      map(seq_len(nsims), function(sim) {
        EpiEstim:::estimate_joint(
          incid[sim, , ,], si_est, priors, seed = 1,
          t_min = 2L, t_max = as.integer(tmax),
          mcmc_control = mcmc_controls
        )
      }
      )
    }
    )
  }
)

