## For 1 location now, edit for 2 locations later.
initial_incidence <- function(type = "growing", n_loc = 1L) {
  ref_i0 <- ifelse(type == "growing", 20, 2000)
  if (n_loc < 2) {
    out <- list(
      incidence::incidence(rep(seq(1, 10), each = ref_i0)),
      incidence::incidence(dates = 10, first_date = 1, last_date = 10)
    )
  } else {
    ## Assume 2 locations, refactor later if needed.
    out <- list(
      incidence::incidence(rep(seq(1, 10), each = ref_i0)),
      incidence::incidence(dates = 10, first_date = 1, last_date = 10),
      ## Location 2, wildtype
      incidence::incidence(rep(seq(1, 10), each = ref_i0)),
      ## Location 2, variant
      incidence::incidence(dates = 10, first_date = 1, last_date = 10)
    )
  }
  out
}
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

  init_days <- nrow(incid_init[[1]]$counts)
  incid <- array(
    NA, dim = c(ndays + init_days, nlocations, nvariants)
  )
  for (loc in seq_len(nlocations)) {
    for (v in seq_len(nvariants)) {
      incid[, loc, v] <- rbind(
        incid_init[[v]]$counts,
        as.matrix( #
          project(
            incid_init[[v]],
            ## R in the future so removing time of seeding
            R = rmatrix[-1, loc, v],
            si = simatrix[, v],
            n_sim = 1,
            n_days = ndays,
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
  ## having this as 20 and starting with 1 case of the variant can lead to an infinite loop
  min_var_cases <- 5
  ## Calculate reproduction number for variant
  rt_variant <- epsilon * rt_ref
  ## Assume reproduction number remains the same
  ## over the time period
  ## Make a vector that goes across rows
  R <- array(NA, dim = c(ndays, n_loc, n_v))
  R[,,1] <- rep(rt_ref, each = ndays)
  R[,,2] <- rep(rt_variant, each = ndays)
  out <- rerun(
    nsims,
    simulate_incidence(
      incid_init, n_loc, n_v, ndays, R, si
    )
  )
  ## total number of cases at the end of the first 10 days
  ## simulation for the variant.
  ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
  out <- out[ncases > min_var_cases]
  message("# of simulations with more than 5 cases ", length(out))
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
    out <- out[ncases > min_var_cases]
    success <- length(out)
  }
  names(out) <- seq_len(nsims)
  out
}

## Same as above but for 2 locations. Since we only
## simulate for 1 or 2 locations, maybe too much of
## an effort to write a wrapper to deal with any
## number of locations but if we create
## simulate_incid_wrapper3, we should refactor the
## function!
## rt_ref2 is the reproduction number for reference in
## location 2.
simulate_incid_wrapper2 <- function(rt_ref, epsilon, si, incid_init,
                                    n_loc = 2, n_v = 2,
                                    ndays = 100, nsims = 100, rt_ref2 = NULL) {
  ## having this as 20 and starting with 1 case of the variant can lead to an infinite loop
  min_var_cases <- 5
  ## Calculate reproduction number for variant
  rt_variant1 <- epsilon * rt_ref
  if (is.null(rt_ref2)) rt_variant2 <- rt_variant1
  else rt_variant2 <- epsilon * rt_ref2
  ## Assume reproduction number remains the same
  ## over the time period
  ## Make a vector that goes across rows
  R <- array(NA, dim = c(ndays, n_loc, n_v))
  R[, 1, 1] <- rep(rt_ref, each = ndays)
  ## For similar locations scenario, rt_ref2 is NULL.
  if (is.null(rt_ref2)) R[, 2, 1] <- rep(rt_ref, each = ndays)
  else R[, 2, 1] <- rep(rt_ref2, each = ndays)
  R[, 1, 2] <- rep(rt_variant1, each = ndays)
  R[, 2, 2] <- rep(rt_variant2, each = ndays)
  out <- rerun(
    nsims,
    simulate_incidence(
      incid_init, n_loc, n_v, ndays, R, si
    )
  )
  ## total number of cases at the end of the first 10 days
  ## simulation for the variant.
  ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
  out <- out[ncases > min_var_cases]
  message("# of simulations with more than 5 cases ", length(out))
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
    out <- out[ncases > min_var_cases]
    success <- length(out)
  }
  names(out) <- seq_len(nsims)
  out
}


## Wrapper for simulating stepwise incidence
simulate_stepwise_incid_wrapper <- function(rt_ref, rt_post_step, step_time,
                                            epsilon, si, incid_init,
                                            n_loc = 1, n_v = 2,
                                            ndays = 100, nsims = 100) {
  ## having this as 20 and starting with 1 case of the variant can lead to an infinite loop
  min_var_cases <- 5
  ## Calculate reproduction number for variant
  rt_variant <- epsilon * rt_ref
  rt_variant_post_step <- epsilon * rt_post_step
  ## Assume reproduction number changes after t = step_time
  ## Make a vector that goes across rows
  R <- array(NA, dim = c(ndays, n_loc, n_v))
  R[,,1] <- c(rep(rt_ref, each = step_time),
              rep(rt_post_step, each = ndays - step_time))
  R[,,2] <- c(rep(rt_variant, each = step_time),
              rep(rt_variant_post_step, each = ndays - step_time))
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
  out <- out[ncases > min_var_cases]
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
    ## check if new simulations have enough cases
    ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
    out <- out[ncases > min_var_cases]
    success <- length(out)
  }
  names(out) <- seq_len(nsims)
  out
}


## Wrapper for simulating stepwise incidence in two locations
simulate_stepwise_incid_wrapper2 <- function(rt_ref, rt_post_step, step_time_l1, step_time_l2,
                                            epsilon, si, incid_init,
                                            n_loc = 2, n_v = 2,
                                            ndays = 100, nsims = 100) {
  ## having this as 20 and starting with 1 case of the variant can lead to an infinite loop
  min_var_cases <- 5
  ## Calculate reproduction number for variant
  rt_variant <- epsilon * rt_ref
  rt_variant_post_step <- epsilon * rt_post_step
  ## Assume reproduction number changes after t = step_time
  ## Make a vector that goes across rows
  R <- array(NA, dim = c(ndays, n_loc, n_v))
  R[,1,1] <- c(rep(rt_ref, each = step_time_l1),
              rep(rt_post_step, each = ndays - step_time_l1))
  R[,2,1] <- c(rep(rt_ref, each = step_time_l2),
               rep(rt_post_step, each = ndays - step_time_l2))
  R[,1,2] <- c(rep(rt_variant, each = step_time_l1),
              rep(rt_variant_post_step, each = ndays - step_time_l1))
  R[,2,2] <- c(rep(rt_variant, each = step_time_l2),
               rep(rt_variant_post_step, each = ndays - step_time_l2))
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
  out <- out[ncases > min_var_cases]
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
    ## check if new simulations have enough cases
    ncases <- map_dbl(out, function(x) sum(x[1:20, , 2]))
    out <- out[ncases > min_var_cases]
    success <- length(out)
  }
  names(out) <- seq_len(nsims)
  out
}
