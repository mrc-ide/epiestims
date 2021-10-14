## orderly::orderly_develop_start(use_draft = "newer", parameters = list(short_run = FALSE))
dir.create("outputs")
simulated_incid <- readRDS("incid.rds")
si_for_est <- readRDS("si_for_est.rds")
param_grid <- readRDS("param_grid.rds")

## Keep subset of simulated_incid
## Epsilon = 2
## Rt = 1.6
## si_mu_variant = 5.4 (i.e multiplier = 1)

incid_subset <- which(param_grid$rt_ref == 1.6 &
                      param_grid$epsilon == 2 &
                      param_grid$si_mu_variant == 5.4)

simulated_incid <- simulated_incid[incid_subset]

si_for_est <- si_for_est[incid_subset]


## Set up epsilon priors to explore

priors <- EpiEstim:::default_priors()
## Change the default priors in EpiEstim
new_default <- epitrix::gamma_mucv2shapescale(1, 1)

priors$epsilon$shape <- new_default$shape
priors$epsilon$scale <- new_default$scale

## Set the shifted epsilon prior
shifted_eps_prior <- epitrix::gamma_mucv2shapescale(1.5, 0.1)
shifted_priors <- priors
shifted_priors$epsilon$shape <- shifted_eps_prior$shape
shifted_priors$epsilon$scale <- shifted_eps_prior$scale

## Combine both sets of priors in list
## first element is default, second is informative

vary_priors <- list(priors, shifted_priors)


mcmc_controls <- list(
  n_iter = 15000L, burnin = 7500L, thin = 20L
  # n_iter = 1500L, burnin = 750L, thin = 10L # test controls
)

tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all
max_attempts <- 3
## Estimate epsilon
plan(multicore)
pwalk(
  list(
    incid = simulated_incid, si = si_for_est,
    index = seq_along(simulated_incid),
    priors = vary_priors, p_index = seq_along(vary_priors)
  ),
  function(incid, si, index, priors, p_index) {
    res <- map(
      tmax_all, function(tmax) {
        ## Loop over the first dimension which is
        ## the set of simulations
        message("tmax = ", tmax)
        future_imap(incid, function(x, i) {
          message("sim = ", i)
          t_min <- EpiEstim::compute_t_min(x, si)
          t_max <- as.integer(t_min + tmax)
          t_max <- min(t_max, nrow(x))
          out <- estimate_joint(
            x, si, priors = priors, seed = 1,
            t_min = t_min,
            t_max = t_max,
            mcmc_control = mcmc_controls
          )
          attempt <- 1
          ## if convergence is achieved, out[["convergence"]] is TRUE
          while (! out[["convergence"]]) {
            message("Attempt ", attempt)
            message("Not yet converged")
            mcmc_controls <- lapply(
              mcmc_controls, function(x) x * 2L
            )
            out <- estimate_joint(
              x, si, priors, seed = 1,
              t_min = t_min,
              t_max = t_max,
              mcmc_control = mcmc_controls
            )
            attempt <- attempt + 1
            ## so that we don't end in an infinite loop
            if (attempt > max_attempts) {
              message("Aborting after 3 attempts")
              ## return whatever you've got.
              return(list(out, out[["convergence"]]))
            }
          }
          list(out, out[["convergence"]])
        }, .options = furrr_options(seed = TRUE, stdout = TRUE),
        .progress = TRUE
        )
      }
    )
    saveRDS(
      res, glue("outputs/estimate_joint_{index}_priors_{p_index}.rds")
    )
  }
)

files2zip <- dir('outputs', full.names = TRUE)
zip(zipfile = "estimate_joint_output.zip", files = files2zip)

