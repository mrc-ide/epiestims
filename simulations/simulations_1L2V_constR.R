# Scenario 1:
# 1 loc, 2 var, constant R (either 1.2 or 3)
require(EpiEstim)
require(incidence)
require(projections)
require(purrr)
require(dplyr)
require(ggplot2)
source("simulation_functions.R")

set.seed(42)

ndays <- 100
n_loc <- 1
n_v <- 2

# SI distr
si_mean <- 6.83
si_std <- 3.8
si <- discr_si(0:30, mu = si_mean, sigma = si_std)
si <- si / sum(si)
si_no_zero <- si[-1]

# Number of simulations
short_run <- TRUE
nsims <- ifelse(short_run, 10, 100)

## Other common things
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 5000L, burnin = as.integer(floor(5e3 / 2)), # speed up
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
  rt_ref = 3,  #c(1.2, 3),
  epsilon = c(1.2, 1.6, 2) #c(seq(from = 1, to = 2, by = 0.1), 2.5, 3)
  #tmax = seq(30, 60, by = 30) # to speed things up
)

index <- seq_len(nrow(sim_params))
sim_params <- sim_params[index, ]

## Remember to use the SI without the first element
## for simulating data (si_no_zero), and to use SI with the
## first element (si) in calls to EpiEstim

si_for_sim <- cbind(si_no_zero, si_no_zero)

si_for_est <- cbind(si, si)

##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
simulated_incid <- pmap(
  list(
    rt_ref = sim_params$rt_ref,
    epsilon = sim_params$epsilon
  ),
  function(rt_ref, epsilon) {
    ## Calculate reproduction number for variant
    rt_variant <- epsilon * rt_ref
    ## Assume reproduction number remains the same
    ## over the time period
    ## Make a vector that goes across rows
    R <- array(NA, dim = c(ndays, n_loc, n_v))
    R[,,1] <- rep(rt_ref, each = ndays)
    R[,,2] <- rep(rt_variant, each = ndays)
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

#saveRDS(simulated_incid, "results/1L2V_incid1.rds") # Rtref 1.2, eps 1, 2 sims 
# saveRDS(simulated_incid, "results/1L2V_incid2.rds") # Rtref 3, eps 2, 10 sims 

tmax_all <- c(30, 40, 50, 60)
names(tmax_all) <- tmax_all

## Estimate epsilon
results <- pmap(
  list(
    incid = simulated_incid
    # tmax = 30
  ),
  function(incid, si) {
    
    map(tmax_all, function(tmax) {
    message("tmax = ", tmax)
    ## Loop over the first dimension which is
    ## the set of simulations
    map(incid, function(x) {
      EpiEstim:::estimate_joint(
        x, si_for_est, priors, seed = 1,
        t_min = 2L, t_max = as.integer(tmax),
        mcmc_control = mcmc_controls
      )
    }
    )
  })
  }
)

params <- as.list(sim_params)
params <- append(
  x = params, values = list(result = results)
)

summary_epsilon <- map_depth(
  results, 3, summarise_epsilon
)
summary_epsilon <- map_depth(
  summary_epsilon, 2, ~ bind_rows(., .id = "sim")
)

summary_epsilon <- map(
  summary_epsilon, ~ bind_rows(., .id = "tmax")
)

summary_epsilon <- map(
  summary_epsilon, function(x) {
    x <- x %>%
      mutate(tmax = recode(tmax,
                           "1" = "30",
                           "2" = "60"))
    x
  }
)


params <- as.list(sim_params)
params <- append(
  params, list(est_epsilon = summary_epsilon)
)

out <- pmap_dfr(
  params, function(rt_ref, epsilon, tmax, est_epsilon) {
    x <- data.frame(
      rt_ref = rt_ref, true_epsilon = epsilon#,
      # tmax = 30
    )
    cbind(x, est_epsilon)
  }
)


if (! dir.exists("results")) dir.create("results")
saveRDS(results, "results/1L2V_raw.rds")
saveRDS(out, "results/1L2V_processed.rds")


### Figures and summary
## out <- readRDS("results/vary_si_processed.rds")

## By rt_ref
ggplot(out) +
  # ylim(1,2.4) +
  geom_point(
    aes(true_epsilon, `50%`)#, position = "dodge2"
  ) +
  geom_linerange(
    aes(true_epsilon, ymin = `2.5%`, ymax = `97.5%`)#,
    # position = "dodge2"
  ) +
  facet_grid(rt_ref~tmax, scales = "free") +
  theme_minimal()
