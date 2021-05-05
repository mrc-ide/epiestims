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
##tmax_all <- seq(ndays, 40, -20)
tmax_all <- 100
tmax_all <- as.integer(tmax_all)
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
  n_iter = 500000, burnin = as.integer(floor(1e4 / 2)),
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
## Diagnostics
## Code to check projections with plot (needs updating with correct variable names)
## TO DO: generate summary grid of incidence plots for each transmission advantage explored

iwalk(
  simulated_incid, function(incid, epsilon) {
    pdf(glue::glue("figures/incid_{epsilon}.pdf"))
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
    dev.off()

  }
)

## Run simulations across all combinations of transmission_advantage and tmax_all
## TO DO: convert code below into function that takes various variables above as arguments
## then we can also look at changes in si
results <- map(
  simulated_incid, function(incid) {
  ## How does varying amount of data affect estimates?
    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)
      EpiEstim:::estimate_joint(
       incid, si_est, priors, seed = 1,
       t_min = 2L, t_max = tmax,
       mcmc_control = mcmc_controls
       )
    }
  )
 }
)
## Weirdly estimate seems to be poorer when tmax is large
## could be convergence issue??
r200 <- results[[1]][[1]]
e200 <- r200$epsilon
rt <- r200[["R"]]

r40 <- results[[1]][[length(tmax_all)]]
e40 <- r40$epsilon
x <- data.frame(tmax = "tmax = 100", eps = e200)
y <- data.frame(tmax = "tmax = 40", eps = e40)
z <- rbind(x, y)

p <- ggplot(x) +
  geom_line(aes(1:nrow(x), eps, col = tmax)) +
  facet_wrap(~tmax, nrow = 2, scales = "free_x") +
   theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

murt <- apply(rt, c(2, 3), mean, na.rm = TRUE)

ggplot() +
  geom_line(aes(1:ncol(murt), murt[2, ])) +
   theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

cowplot::save_plot("figures/possible_convergence_issue.pdf", p)
## Maximum
## cumulative incidence at tmax across
## locations variants
cum_incid <- map_dfr(
  simulated_incid, function(incid) {
    map_dfr(tmax_all, function(tmax) {
      x <- incid[seq_len(tmax), ,]
      out <- apply(x, c(2, 3), sum)
      data.frame(
        max_cum_incid = max(out)
      )
    }, .id = "tmax"
    )
  }, .id = "epsilon"
)

vary_tmax_est <- map_depth(
  results, 2, process_fit
)

vary_tmax_est <- map_dfr(
  vary_tmax_est, function(x) {
    bind_rows(x, .id = "tmax")
  }, .id = "epsilon"
)

vary_tmax_est <- left_join(vary_tmax_est, cum_incid)

vary_tmax_est$tmax <- factor(
  vary_tmax_est$tmax, levels = rev(tmax_all), ordered = TRUE
)
## Duplicate column so that we can plot true values
vary_tmax_est$true_epsilon <- as.numeric(vary_tmax_est$epsilon)

vary_tmax_est$epsilon <- factor(
  vary_tmax_est$epsilon,
  levels = transmission_advantage,
  ordered = TRUE
)


est_epsilon <- vary_tmax_est[vary_tmax_est$param == "epsilon", ]


p1 <- ggplot(est_epsilon) +
  geom_linerange(
    aes(tmax, ymin = mu - sd, ymax = mu + sd)
  ) +
  geom_point(aes(tmax, mu)) +
  geom_hline(
    aes(yintercept = true_epsilon),
    col = "red", linetype = "dashed"
  ) +
  expand_limits(y = 0) +
  ylab("epsilon") +
  facet_wrap(~ epsilon, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  )

cowplot::save_plot("figures/epsilon_tmax.pdf", p1)

p2 <- ggplot(est_epsilon) +
  geom_linerange(
    aes(log(max_cum_incid), ymin = mu - sd, ymax = mu + sd)
  ) +
  geom_point(aes(log(max_cum_incid), mu)) +
  geom_hline(
    aes(yintercept = true_epsilon),
    col = "red", linetype = "dashed"
  ) +
  facet_wrap(~ epsilon, ncol = 2, scales = "free_x") +
  expand_limits(y = 0) +
  ylab("epsilon") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  )

cowplot::save_plot("figures/epsilon_cumincid.pdf", p2)

