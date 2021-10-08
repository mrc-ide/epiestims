## orderly::orderly_develop_start(use_draft = "newer", parameters = list(short_run = FALSE))
dir.create("outputs")
simulated_incid <- readRDS("incid.rds")
si_for_est <- readRDS("si_for_est.rds")

priors <- EpiEstim:::default_priors()
mcmc_controls <- list(
  n_iter = 15000L, burnin = 7500L, thin = 20L
)

tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all
max_attempts <- 3
## Estimate epsilon
plan(multicore)
iwalk(
  simulated_incid,
  function(incid, index) {
    res <- map(
      tmax_all, function(tmax) {
        ## Loop over the first dimension which is
        ## the set of simulations
        message("tmax = ", tmax)
        future_imap(incid, function(x, i) {
          message("sim = ", i)
          t_min <- EpiEstim::compute_t_min(x, si_for_est)
          t_max <- as.integer(t_min + tmax)
          t_max <- min(t_max, nrow(x))
          out <- estimate_joint(
            x, si_for_est, priors, seed = 1,
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
              x, si_for_est, priors, seed = 1,
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
        }, .options = furrr_options(seed = TRUE, stdout = FALSE),
        .progress = TRUE
        )
      }
    )
    saveRDS(
      res, glue("outputs/estimate_joint_{index}.rds")
    )
  }
)

files2zip <- dir('outputs', full.names = TRUE)
zip::zip(zipfile = "estimate_joint_output.zip", files = files2zip)

