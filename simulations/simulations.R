## Things we want to test
## R levels: <1, 1.5, 4.5
## X
## Transmission Advantage: low, medium, high say
## 0.1, 0.5, 1

## 1. Try single location first.
## Update - doesn't work.
library(EpiEstim)
library(projections)
library(incidence)
library(purrr)
library(dplyr)
library(ggplot2)

source("simulations/simulation_functions.R")

seed <- 42
set.seed(seed)
ndays <- 200
n_loc <- 3
n_v <- 2
## Set-up reference reproduction numbers
rt_ref <- c(1.5, 1.2, 1.1)
transmission_advantage <- seq(1, 2, 0.1)

results <- map(transmission_advantage, function(advantage) {
  
  message("transmission advantage = ", advantage)
  
  ## Reproduction number for variant
  rt_variant <- advantage * rt_ref
  ## Assume reproduction number remains the same
  ## over the time period
  ## Make a vector that goes across rows
  R <- array(NA, dim = c(ndays, n_loc, n_v))
  R[,,1] <- rep(rt_ref, each = ndays) 
  R[,,2] <- rep(rt_variant, each = ndays)
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
  ## start with 20 infected individual
  ## large-ish initial seed to ensure we get an epidemic
  ## can explore lower seed numbers later
  initial_incidence <- incidence::incidence(rep(1, 20))
  
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
  
  # check projections with plot
  plot(log(1 + incid[, 1, 1]), type= "l", xlab = "date", ylab = "log(1 + Incidence)")
  lines(log(1 + incid[, 2, 1]), col = "blue")
  lines(log(1 + incid[, 3, 1]), col = "red")
  lines(log(1 + incid[, 1, 2]), lty = 2)
  lines(log(1 + incid[, 2, 2]), col = "blue", lty = 2)
  lines(log(1 + incid[, 3, 2]), col = "red", lty = 2)
  legend("bottomright", c("Strain 1, location 1", "Strain 1, location 2", "Strain 1, location 3",
                          "Strain 2, location 1", "Strain 2, location 2", "Strain 2, location 3"),
         col = c("black", "blue", "red", "black", "blue", "red"),
         lty = c(1, 1, 1, 2, 2, 2), cex = 0.7)
  
  ## Now estimate epsilon
  priors <- EpiEstim:::default_priors()
  mcmc_controls <- list(n_iter = 1e4L, burnin = as.integer(floor(1e4 / 2)), thin = 10L)
  
  ## x is a list with elements R and epsilon
  ## epsilon is  a vector with length (n_iter - burnin)/thin + 1
  ## R is an array with dimensions T X n_loc X (n_iter - burnin)/thin + 1
  
  
  ## How does varying amount of data affect estimates?
  tmax_all <- as.integer(seq(200, 50, -10))
  names(tmax_all) <- tmax_all
  vary_tmax_fits <- map(
    tmax_all, function(tmax) {
      message("tmax = ", tmax)
      EpiEstim:::estimate_joint(
        incid, si_distr, priors, seed = 1,
        t_min = 2L, t_max = tmax,
        mcmc_control = mcmc_controls
      )
    }
  )
  
  vary_tmax_est <- imap_dfr(
    vary_tmax_fits, function(out, tmax) {
      message("tmax = ", tmax)
      process_fit(out)
    }, .id = "tmax"
  )
  
  vary_tmax_est$tmax <- factor(
    vary_tmax_est$tmax, levels = rev(tmax_all), ordered = TRUE
  )
  
  vary_tmax_est
  
})


names(results) <- transmission_advantage
results_comb <- dplyr::bind_rows(results, .id = "transmission_advantage")
results_comb$transmission_advantage <- as.numeric(results_comb$transmission_advantage)

trans_adv <- results_comb %>% 
  group_by(transmission_advantage) %>% 
  summarise(adv = mean(transmission_advantage))

ggplot(data = results_comb[results_comb$param == "epsilon", ]) +
  geom_linerange(aes(tmax, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_point(aes(tmax, `50%`)) +
  geom_hline(data = trans_adv, aes(
    yintercept = adv),
    linetype = "dashed", color = "red"
  ) +
  facet_wrap(~ transmission_advantage, nrow = 3, ncol = 4)

