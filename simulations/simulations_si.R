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
## Number of simulations for each parameter
nsims <- 10

## Reference reproduction numbers in each location
##rt_ref <- c(1.5, 1.5, 1.5)
rt_ref <- c(1.5, 1.5)

## Range of transmission advantage values to explore
## transmission_advantage <- seq(2, 3, 0.2)
transmission_advantage <- 1.2
names(transmission_advantage) <- transmission_advantage

## Define range of tmax values to explore
## tmax_all <- seq(ndays, 40, -20)
tmax_all <- c(ndays, 40)
tmax_all <- as.integer(tmax_all)
names(tmax_all) <- tmax_all

## Set serial interval distributions
## TO DO: vary these (eg variant has si that is shorter or longer than ref)
## covid serial interval from IBM
##si_mean <- c(3.41, 6.83, 13.66)
si_mean_ref <- 6.83
## Variant SI
si_multiplier <- c(0.25, 0.5, 1, 1.5, 2)
si_mean <- si_mean_ref * si_multiplier
names(si_mean) <- si_mean
si_std <- 3.8

si_variant <- map(si_mean, function(x) {
  ##shape_scale <- epitrix::gamma_mucv2shapescale(x, si_std / x)
  ##cutoff <- qgamma(0.99, shape_scale$shape, scale = shape_scale$scale)
  ## 30 seems to capture most of the probability mass alright
  si <- discr_si(0:30, mu = x, sigma = si_std)
  si <- si / sum(si)
})
ref_index <- match(si_mean_ref, names(si_mean))
si_ref <- si_variant[[ref_index]] # mean 6.83 as per covid IBM
si_no_zero_ref <- si_ref[-1]
si_no_zero_var <- map(si_variant, function(x) x[-1])
## Needs to change if number of variants is
## changed.
## For simulation
si_distr <- map(si_no_zero_var, function(x) cbind(si_no_zero_ref, x))
## For estimation
si_est <- map(
  si_variant, function(x) cbind(si_ref, x)
)
priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 10000L, burnin = as.integer(floor(1e4 / 2)),
  thin = 10L
)


simulated_incid <- imap(
  si_distr, function(x, si_name) {
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
      initial_incidence, n_loc, n_v, ndays, R, x, nsims = nsims
    )
  })

## Vary SI of variant (2xsmaller, 2xlarger and same as ref) for variant with transm advantage of 2
## In estimating epsilon here, we assume that we know the SI
results <- imap(
  simulated_incid, function(incid, si_name) {
    message("si_mean_var = ", si_name)
    map(tmax_all, function(tmax) {
      message("tmax = ", tmax)
      ## Loop over the first dimension which is
      ## the set of simulations
      map(seq_len(nsims), function(sim) {
        EpiEstim:::estimate_joint(
           incid[sim, , ,], si_est[[si_name]], priors, seed = 1,
           t_min = 2L, t_max = as.integer(tmax),
           mcmc_control = mcmc_controls
        )
      }
      )
    }
    )
  }
)
if (! dir.exists("results")) dir.create("results")
saveRDS(results, "results/vary_si.rds")

## results are now in the 3rd level.
## That is, results for incidence simulation 1
## for si 1, tmax 1, are in results[[1]][[1]]
## results for incidence simulation 2
## for si 1, tmax 1, are in results[[1]][[2]]
vary_si <- map_depth(results, 3, process_fit)

out <- map_depth(vary_si, 2, function(x) imap_dfr(x, ~ .x, .id = "sim"))
out2 <- map(out, function(x) map_dfr(x, ~ .x, .id = "tmax"))
vary_si <- map_dfr(out2, function(x) x, .id = "si")

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
  geom_point(aes(si, mu), alpha = 0.3) +
  geom_hline(
    aes(yintercept = true_epsilon),
    col = "red", linetype = "dashed", alpha = 0.3
  ) +
  expand_limits(y = 0) +
  scale_x_discrete(
    breaks = levels(vary_si$si),
    labels = glue::glue(" X {si_multiplier}")
  ) +
  xlab("Variant serial interval") +
  ylab("epsilon") +
  facet_wrap(~ tmax, ncol = 2, labeller = label_both) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  )

p1

cowplot::save_plot("figures/epsilon_tmax_si.pdf", p1)

## summary
x <- group_by(est_epsilon, si, tmax) %>%
  summarise(
    in_95CrI = sum(transmission_advantage > `2.5%` & transmission_advantage < `97.5%`)
  )


