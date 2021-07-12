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
