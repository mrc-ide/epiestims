## Things we want to test
## R levels: <1, 1.5, 4.5
## X
## Transmission Advantage: low, medium, high say
## 0.1, 0.5, 1

## 1. Try single location first
library(EpiEstim)
library(projections)
library(incidence)
seed <- 42
set.seed(seed)
ndays <- 100
n_loc <- 1
n_v <- 2
R1_loc1 <- 1.5
transmission_advantage <- 1.1
R2_loc1 <- transmission_advantage * R1_loc1


R <- array(NA, dim = c(ndays, n_loc, n_v))
R[, 1, 1] <- R1_loc1
R[, 1, 2] <- R2_loc1

#####################################################################
### Simulate epidemics with those reproduction numbers ###
#####################################################################
## covid serial interval from IBM
si_mean <- 6.83
si_std <- 3.8
si <- EpiEstim::discr_si(0:30, mu = si_mean, sigma = si_std)
si <- si / sum(si)
si_no_zero <- si[-1]
si_distr <- cbind(si, si)
## start with 5 infected individual
initial_incidence <- incidence::incidence(rep(1, 5))

incid <- array(NA, dim = c(ndays, n_loc, n_v))

for (loc in 1:n_loc) {
  for (v in 1:n_v) {
    incid[, loc, v] <- rbind(
      initial_incidence$counts,
      as.matrix( #
        projections::project(initial_incidence,
          R = R[-1, loc, v], # R in the future so removing time of seeding
          si = si_no_zero,
          n_sim = 1,
          n_days = ndays - 1,
          time_change = seq_len(length(R[, loc, v]) - 2) - 1
        )
      )
    )
  }
}

## Now estimate epsilon
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(n_iter = 1e4L, burnin = floor(1e4 / 2), thin = 10L)
x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1, t_min = 2L, t_max = 100L, mcmc_control = mcmc_controls)
