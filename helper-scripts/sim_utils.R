##' Simulate incidence for multiple locations and multiple
##' variants
##' No checks implemented, make sure you input right things in
##' right dimensions
##' @param incid_init initial incidence as a list of incidence
##' objects. Each list element is the initial incidence object
##' for a variant.
##' @param nlocations number of locations
##' @param nvariants number of variants
##' @param ndays number of days
##' @param rmatrix matrix of reproduction numbers to use
##' dimension: ndays X nlocations X nvariants
##' @param simatrix matrix of SI distribution, 1 column for
##' each variant
##' @param nsims number of simulations, defaults to 1
##' @return matrix of
##' dimensions nsims X ndays X nlocations X nvariants
##' @author Sangeeta Bhatia, Jack Wardle
simulate_incidence <- function(incid_init, nlocations,
                               nvariants, ndays, rmatrix,
                               simatrix) {
  incid <- array(NA, dim = c(ndays, nlocations, nvariants))

  for (loc in seq_len(nlocations)) {
    for (v in seq_len(nvariants)) {
      incid[, loc, v] <- rbind(
        tail(incid_init[[v]]$counts, 1),
        as.matrix( #
          project(
            incid_init[[v]],
            ## R in the future so removing time of seeding
            R = rmatrix[-1, loc, v],
            si = simatrix[, v],
            n_sim = 1,
            n_days = ndays - 1,
            instantaneous_R = TRUE,
            time_change = seq_len(
              length(rmatrix[, loc, v]) - 2
            ) - 1
          )
        )
      )
    }
  }
  incid
}


simulate_incid_wrapper <- function(rt_ref, epsilon, si, incid_init,
                                   n_loc = 1, n_v = 2,
                                   ndays = 100, nsims = 100) {
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
    out <- rerun(
      nsims,
      simulate_incidence(
        incid_init, n_loc, n_v, ndays, R, si
      )
    )
  ## total number of cases at the end of the first 10 days
  ## simulation for the variant.
  ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
  out <- out[ncases > 20]
  message("# of simulations with more than 20 cases ", length(out))
  ## At this point out may have less than
  ## the desired number of simulations
  success <- length(out)
  while (success < nsims) {
    message("More sims needed ", nsims - success)
    more <- rerun(
      nsims - success,
      simulate_incidence(
        incid_init, n_loc, n_v, ndays, R, si
      )
    )
    out <- append(out, more)
    ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
    out <- out[ncases > 20]
    success <- length(out)
  }
  names(out) <- seq_len(nsims)
  out
}

## nicked from EpiEstim@fix_tmin; After PR 127 is
## merged, this function can be deleted.
compute_si_cutoff <- function(si_distr, miss_at_most = 0.05) {
  if (sum(si_distr) != 1) {
    warning("Input SI distribution should sum to 1. Normalising now")
    si_distr <- si_distr / colSums(si_distr)
  }
  cutoff <- 1 - miss_at_most
  cdf <- apply(si_distr, 2, cumsum)
  idx <- apply(
    cdf, 2,
    function(col) Position(function(x) x > cutoff, col)
  )
  as.integer(max(idx))
}
