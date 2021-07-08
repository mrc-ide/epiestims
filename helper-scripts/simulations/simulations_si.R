## Serial interval mean sensitivity as 1 location, stepwise R
## But with SI for variant 2 which has mean x
## times the mean SI for reference variant with
## x = c(0.5, 0.75, 1, 1.25, 1.5)
## Estimation is done with
## 1. Correct SI
## 2. Incorrect SI (of reference variant)
## R_reference c(1.2, 3)
## epsilon c(seq(from = 1, to = 2, by = 0.1), 2.5, 3)
if (! dir.exists("results")) dir.create("results")
source("global.R")
short_run <- FALSE
ndays <- 100
n_loc <- 1
n_v <- 2
## Number of data sets simulated for each parameter

si_mu_ref <- 6.83
si_std_ref <- 3.8

## Other common things
priors <- list(epsilon = list(shape = 1, scale = 1), 
               R = list(shape = 1, scale = 1))
mcmc_controls <- list(
  n_iter = 10000L, burnin = as.integer(floor(1e4 / 2)),
  thin = 10L
)

## Seed with 1 case but select simulations
## that went on to generate at least 20 cases
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

tmax_all <- seq(20, 60, 10)
names(tmax_all) <- tmax_all

index <- if (short_run) c(1, nrow(sim_params)) else seq_len(nrow(sim_params))
nsims <- ifelse(short_run, 10, 100)
sim_params <- sim_params[index, ]


