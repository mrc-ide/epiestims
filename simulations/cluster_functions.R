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
      seq_len(10 * nsims), function(x) {
        simulate_incidence(
          initial_incidence, n_loc, n_v, ndays, R, si
        )
      }
    )
    ## total number of cases at the end of the first 10 days
    ## simulation for the variant.
    ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
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
      out[[index]][[sim]] <-
        EpiEstim::estimate_joint(
                 incid[[sim]], si_for_est, priors, seed = 1,
                 t_min = 2L, t_max = as.integer(tmax),
                 mcmc_control = mcmc_controls
                 )

    }
  }
  out
}


manager <- function(rt_ref, epsilon, si_mu_variant, si_std_variant) {
  si_distr_ref <- EpiEstim::discr_si(
    0:30, mu = si_mu_ref, sigma = si_std_ref
  )
  si_no_zero_ref <- si_distr_ref[-1]
  si_distr_variant <- EpiEstim::discr_si(
    0:30, mu = si_mu_variant, sigma = si_std_variant
  )
  si_distr_variant <- si_distr_variant / sum(si_distr_variant)
  si_no_zero_var <- si_distr_variant[-1]
  si_for_sim <- cbind(si_no_zero_ref, si_no_zero_var)
  ## si_for_est <- cbind(si_distr_ref, si_distr_variant)
  ## Mis-specify SI for variant
  si_for_est <- cbind(si_distr_ref, si_distr_ref)
  incid <- simulate_incid_wrapper(rt_ref, epsilon, si_for_sim)
  res <- estimate_wrapper(incid, si_for_est)
  list(incid = incid, res = res)
}
