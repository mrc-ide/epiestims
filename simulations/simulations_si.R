## Things we want to test
## R levels: <1, 1.5, 4.5
## X
## Transmission Advantage: low, medium, high say
## 0.1, 0.5, 1

## 1. Try single location first.
## Update - doesn't work.

## Common parameters
## Set time, locations and number of variants in simulations
source("global.R")
ndays <- 100
n_loc <- 2
n_v <- 2

## Reference reproduction numbers in each location
##rt_ref <- c(1.5, 1.5, 1.5)
rt_ref <- c(1.5, 1.5)

## Range of transmission advantage values to explore
## transmission_advantage <- seq(2, 3, 0.2)
transmission_advantage <- 1.2
names(transmission_advantage) <- transmission_advantage

## Define range of tmax values to explore
tmax_all <- seq(ndays, 40, -20)
# tmax_all <- c(200, 40)
tmax_all <- as.integer(tmax_all)
names(tmax_all) <- tmax_all

## Set serial interval distributions
## TO DO: vary these (eg variant has si that is shorter or longer than ref)
## covid serial interval from IBM
##si_mean <- c(3.41, 6.83, 13.66)
si_mean_ref <- 6.83
## Variant SI
si_mean <- si_mean_ref * seq(0.2, 3, 0.2)
names(si_mean) <- si_mean
si_std <- 3.8
si_variant <- map(si_mean, function(x) {
  ##shape_scale <- epitrix::gamma_mucv2shapescale(x, si_std / x)
  ##cutoff <- qgamma(0.99, shape_scale$shape, scale = shape_scale$scale)
  ## 30 seems to capture most of the probability mass alright
  si <- discr_si(0:30, mu = x, sigma = si_std)
  si <- si / sum(si)
})

si_ref <- si_variant[["6.83"]] # mean 6.83 as per covid IBM
si_no_zero_ref <- si_ref[-1]
si_no_zero_var <- map(si_variant, function(x) x[-1])
## Needs to change if number of variants is
## changed.
## For simulation
si_distr <- map(si_no_zero_var, function(x) {
  cbind(si_no_zero_ref, x)
})
## For estimation
si_est <- map(si_variant, function(x) {
  cbind(si_ref, x)
})
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 10000L, burnin = as.integer(floor(1e4 / 2)),
  thin = 10L
)

## TODo
## Simulate incidence for reference variant, so that it
## stays the same across the epsilon values we explore


simulated_incid <- imap(
  si_distr, function(si_distr, si_name) {
    message("variant si_mean = ", si_name)
    ## Calculate reproduction number for variant
    rt_variant <- transmission_advantage * rt_ref
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
## Diagnostics
## Code to check projections with plot (needs updating with correct variable names)
## TO DO: generate summary grid of incidence plots for each transmission advantage explored

iwalk(
  simulated_incid, function(incid, si) {
    plot(log(1 + incid[, 1, 1]), type= "l", xlab = "date", ylab = "log(1 + Incidence)")
    lines(log(1 + incid[, 1, 2]), lty = 2)
    legend("bottomright", c("Strain 1", "Strain 2"),
           lty = c(1, 2), cex = 0.7)
    title(paste0("si_mean_ref = 6.83; si_mean_var = ", si))

  }
)

## Vary SI of variant (2xsmaller, 2xlarger and same as ref) for variant with transm advantage of 2
## In estimating epsilon here, we assume that we know the SI
results <- imap(
  simulated_incid, function(incid, si_name) {

    message("si_mean_var = ", si_name)

    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)

      EpiEstim:::estimate_joint(
        incid, si_est[[si_name]], priors, seed = 1,
        t_min = 2L, t_max = as.integer(tmax),
        mcmc_control = mcmc_controls
      )
    }
    )
  }
)




## Plot estimated values
vary_si <- map_depth(results, 2, process_fit)

vary_si <- map_dfr(
  vary_si, function(x) {
    bind_rows(x, .id = "tmax")
  }, .id = "si"
)

vary_si$tmax <- factor(
  vary_si$tmax, levels = rev(tmax_all), ordered = TRUE
)

vary_si$true_epsilon <- as.numeric(transmission_advantage)

vary_si$si <- factor(
  vary_si$si, levels = si_mean, ordered = TRUE
)


est_epsilon <- vary_si[vary_si$param == "epsilon", ]


p1 <- ggplot(est_epsilon) +
  geom_linerange(
    aes(si, ymin = mu - sd, ymax = mu + sd)
  ) +
  geom_point(aes(si, mu)) +
  geom_hline(
    aes(yintercept = true_epsilon),
    col = "red", linetype = "dashed"
  ) +
  expand_limits(y = 0) +
  ylab("epsilon") +
  facet_wrap(~ tmax, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  )
p1

cowplot::save_plot("figures/epsilon_tmax_si.pdf", p1)


## Now consider a situation where we estimate epsilon assuming that
## we do not know the si of new variant so use current si

results_same_si <- imap(
  simulated_incid, function(incid, si_name) {

    message("si_mean_var = ", si_name)

    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)

      EpiEstim:::estimate_joint(
        incid, si_est[["6.83"]], priors, seed = 1,
        t_min = 2L, t_max = as.integer(tmax),
        mcmc_control = mcmc_controls
      )
    }
    )
  }
)




## Plot the estimated values of epsilon

vary_si_sameref <- map_depth(
  results_same_si, 2, process_fit
)

vary_si_sameref <- map_dfr(
  vary_si_sameref, function(x) {
    bind_rows(x, .id = "tmax")
  }, .id = "si"
)

vary_si_sameref$tmax <- factor(
  vary_si_sameref$tmax, levels = rev(tmax_all), ordered = TRUE
)

vary_si_sameref$true_epsilon <- as.numeric(transmission_advantage)

vary_si_sameref$si <- factor(
  vary_si_sameref$si,
  levels = si_mean,
  ordered = TRUE
)


est_epsilon_sameref <- vary_si_sameref[vary_si_sameref$param == "epsilon", ]


p2 <- ggplot(est_epsilon_sameref) +
  geom_linerange(
    aes(si, ymin = mu - sd, ymax = mu + sd)
  ) +
  geom_point(aes(si, mu)) +
  geom_hline(
    aes(yintercept = true_epsilon),
    col = "red", linetype = "dashed"
  ) +
  expand_limits(y = 0) +
  ylab("epsilon") +
  facet_wrap(~ tmax, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  )
p2

cowplot::save_plot("figures/epsilon_tmax_si_sameref.pdf", p2)

