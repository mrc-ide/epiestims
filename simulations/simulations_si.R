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
priors <- EpiEstim:::default_priors()
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

tmax_all <- seq(10, 60, 10)
names(tmax_all) <- tmax_all

index <- if (short_run) c(1, nrow(sim_params)) else seq_len(nrow(sim_params))
nsims <- ifelse(short_run, 10, 100)
sim_params <- sim_params[index, ]
si_distr_ref <- discr_si(
  0:30, mu = si_mu_ref, sigma = si_std_ref
)
## Remember to use the SI without the first element
## for simulating data, and to use SI with the
## first element in calls to EpiEstim
si_no_zero_ref <- si_distr_ref[-1]
## shape_scale <- epitrix::gamma_mucv2shapescale(x, si_std / x)
## cutoff <- qgamma(0.99, shape_scale$shape, scale = shape_scale$scale)
## 30 seems to capture most of the probability mass alright
si_distr_variant <- pmap(
  sim_params[, c("si_mu_variant", "si_std_variant")], function(si_mu_variant, si_std_variant) {
  si <- discr_si(0:30, mu = si_mu_variant, sigma = si_std_variant)
  si / sum(si)
  }
)

si_no_zero_var <- map(si_distr_variant, function(x) x[-1])
si_for_sim <- map(
  si_no_zero_var, function(x) cbind(si_no_zero_ref, x)
)
si_for_est <- map(
  si_distr_variant, function(x) cbind(si_distr_ref, x)
)


##############################################################################
## Simulate epidemic incidence data with input reproduction numbers and si  ##
##############################################################################
simulated_incid <- pmap(
  list(
    rt_ref = sim_params$rt_ref,
    epsilon = sim_params$epsilon,
    si = si_for_sim
  ),
  function(rt_ref, epsilon, si) {
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
      seq_len(2 * nsims), function(x, index) {
        message("Simulation ", index)
        simulate_incidence(
          initial_incidence, n_loc, n_v, ndays, R, si
        )
      }
    )
    ## total number of cases at the end of the
    ## simulation for the variant.
    ncases <- map_dbl(out, function(x) sum(x[, , 2]))
    out <- out[ncases > 20]
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

saveRDS(simulated_incid, "results/vary_si_simulated_data.rds")
## simulated_incid <- readRDS("results/vary_si_simulated_data.rds")
## In estimating epsilon here, we assume that we know the SI


results <- pmap(
  list(incid = simulated_incid, si = si_for_est),
  function(incid, si) {
    map(tmax_all, function(tmax) {
        message("tmax = ", tmax)
        imap(incid, function(x, index) {
          message("Sim # ", index)
          tryCatch(
          {
            out <- EpiEstim:::estimate_joint(
              x, si, priors, seed = 1,
              t_min = 2L, t_max = as.integer(tmax),
              mcmc_control = mcmc_controls
              )
            list(epsilon = out[["epsilon"]][, seq_len(tmax)])
             }, error = function(cond) {
               out <- list()
               class(out) <- "error"
               out
             }
          )
        }
        )
      }
    )
  }
)
saveRDS(results, "results/vary_si_raw.rds")

params <- as.list(sim_params)
params <- append(
  x = params, values = list(result = results)
)

summary_epsilon <- map_depth(
  results, 3, function(x){
    if (inherits(x, "error")) NULL
    else summarise_epsilon(x)
  }
)
## Get rid of nulls now
summary_epsilon <- map(summary_epsilon, function(incid_eps) {
  map_dfr(
    incid_eps, function(tmax_eps) bind_rows(tmax_eps, .id = "sim"),
    .id = "tmax"
  )
}
)


params <- as.list(sim_params)
params <- append(
  params, list(est_epsilon = summary_epsilon)
)

out <- pmap_dfr(
  params, function(rt_ref, epsilon, si_mu_variant,
                   si_std_variant, est_epsilon) {
    x <- data.frame(
      rt_ref = rt_ref, true_epsilon = epsilon,
      si_mu_variant = si_mu_variant,
      si_std_variant = si_std_variant
    )
    cbind(x, est_epsilon)
  }
)




saveRDS(out, "results/vary_si_processed.rds")


