simulate_incid_wrapper <- function(rt_ref, epsilon, si) {
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
    out <- map(
      seq_len(2 * nsims), function(x) {
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


estimate_wrapper <- function(incid, si_for_est) {

  out <- vector(mode = "list", length = length(tmax_all))
  for (index in seq_along(out)) {
    tmax <- tmax_all[[index]]
    for (sim in seq_along(incid)) {
      out[[index]][[sim]] <-  tryCatch({
        res <- EpiEstim:::estimate_joint(
                 incid[[sim]], si_for_est, priors, seed = 1,
                 t_min = 2L, t_max = as.integer(tmax),
                 mcmc_control = mcmc_controls
                 )
        list(epsilon = res[["epsilon"]][, seq_len(tmax)])
      }, error = function(cond) {
        out <- list()
        class(out) <- "error"
        out
      }
      )
    }
  }
  out
}


manager <- function(rt_ref, epsilon, si_for_sim, si_for_est) {
  incid <- simulate_incid_wrapper(rt_ref, epsilon, si)
  estimate_wrapper(incid, si_for_est)
}
