## Things we want to test
## R levels: <1, 1.5, 4.5
## X
## Transmission Advantage: low, medium, high say
## 0.1, 0.5, 1
## Debugging

## 1. Try single location first.
## Update - doesn't work.
library(EpiEstim)
library(projections)
library(incidence)
library(purrr)
library(dplyr)
library(ggplot2)
library(here)
source("simulations/simulation_functions.R")

if (! dir.exists("figures")) dir.create("figures")
seed <- 42
set.seed(seed)
## Common parameters
## Set time, locations and number of variants in simulations
ndays <- 100
n_loc <- 2
n_v <- 2

## Reference reproduction numbers in each location
##rt_ref <- c(1.5, 1.5, 1.5)
rt_ref <- c(1.5, 1.5)

## Range of transmission advantage values to explore
## transmission_advantage <- seq(2, 3, 0.2)
transmission_advantage <- 2
names(transmission_advantage) <- transmission_advantage

## Define range of tmax values to explore
tmax_all <- c(100, 40) # seq(ndays, 40, -20)
##tmax_all <- 100
##tmax_all <- as.integer(tmax_all)
names(tmax_all) <- tmax_all

## Set serial interval distributions
## TO DO: vary these (eg variant has si that is shorter or longer than ref)
## covid serial interval from IBM
si_mean <- 6.83
si_std <- 3.8
si <- discr_si(0:30, mu = si_mean, sigma = si_std)
si <- si / sum(si)
si_no_zero <- si[-1]
## Needs to change if number of variants is
## changed.
## For simulation
si_distr <- cbind(si_no_zero, si_no_zero)
## For estimation
si_est <- cbind(si, si)
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 500000L, burnin = as.integer(floor(1e4 / 2)),
  thin = 100L
)


## TODo
## Simulate incidence for reference variant, so that it
## stays the same across the epsilon values we explore


simulated_incid <- map(
  transmission_advantage, function(advantage) {
    message("transmission advantage = ", advantage)
    ## Calculate reproduction number for variant
    rt_variant <- advantage * rt_ref
    ## Assume reproduction number remains the same
    ## over the time period
    ## Make a vector that goes across rows

    ## TO DO: update this so that it can handle n_v variants
    R <- array(NA, dim = c(ndays, n_loc, n_v))
    R[,,1] <- rep(rt_ref, each = ndays)
    R[,,2] <- rep(rt_variant, each = ndays)
  ##############################################################################
  ## Simulate epidemic incidence data with input reproduction numbers and si  ##
  ##############################################################################
  ## start with 20 infected individuals
  ## large-ish initial seed to ensure we get an epidemic with lower values of R
  ## can explore lower seed numbers later
    initial_incidence <- incidence::incidence(rep(1, 20))
    simulate_incidence(
      initial_incidence, n_loc, n_v, ndays, R, si_distr
    )
  })



## Run simulations across all combinations of transmission_advantage and tmax_all
## TO DO: convert code below into function that takes various variables above as arguments
## then we can also look at changes in si
results <- map(
  simulated_incid, function(incid) {
  ## How does varying amount of data affect estimates?
    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)
      
      # modify incidence array so that the variant with highest incidence at tmax is the reference
      max_incidence <- arrayInd(which.max(incid[tmax,,]),
                                dim(incid[tmax,,]))
      incid_new_ref <- array(NA,
                             dim = dim(incid))
      incid_new_ref[ , , 1] <- incid[ , , max_incidence[2]]
      incid_new_ref[ , , -1] <- incid[, , -max_incidence[2]]
      
      # now call estimate_joint from EpiEstim using the re-ordered incidence data
      estimate_joint(
       incid_new_ref, si_est, priors, seed = 1,
       t_min = 2L, t_max = as.integer(tmax),
       mcmc_control = mcmc_controls
      )
    }
  )
 }
)

vary_tmax <- map_depth(
  results, 2, function(x) {
    rpost <- x$r_posterior_params
    epost <- x$eps_posterior_params
    emucv <- apply(
      epost, 1, function(row) epitrix::gamma_shapescale2mucv(row[1], row[2])
    )
    emucv <- map_dfr(
      emucv, function(x) data.frame(mu = x$mu, cv = x$cv, sd = x$mu * x$cv)
    )
    emucv$iter <- seq_len(nrow(emucv))
    emucv$low <- qgamma(
      0.025, shape = epost$shape, scale = epost$scale
    )
    emucv$high <- qgamma(
      0.975, shape = epost$shape, scale = epost$scale
    )

    emucv
  }
)

vary_tmax <- map_dfr(
  vary_tmax, function(x) bind_rows(x, .id = "tmax")
)
vary_tmax$tmax <- factor(
  vary_tmax$tmax, levels = tmax_all, ordered = TRUE
)

p <- ggplot(vary_tmax) +
  geom_linerange(
    aes(x = iter, ymin = low, ymax = high)
  ) +
  facet_wrap(~tmax, ncol = 1) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(
  filename = glue::glue("figures/tmax_epsilon_qntls_{rt_ref[1]}_{transmission_advantage}.pdf"),
  p
)
