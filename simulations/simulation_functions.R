summarise_R <- function(fit, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) {
  ## Get rid of the first row because of NAs
  r_est <- apply(
    fit$R[-1, , ], c(1, 2), quantile, probs = probs, na.rm = na.rm
  )
  r_mu <- apply(fit$R[-1, , ], c(1, 2), mean, na.rm = na.rm)
  r_sd <- apply(fit$R[-1, , ], c(1, 2), sd, na.rm = na.rm)
  nt <- dim(r_est)[2]
  nl <- dim(r_est)[3]
  r_estdf <- data.frame(
    time = rep(NA, nt * nl),
    location = rep(NA, nt * nl),
    `2.5%` = rep(NA, nt * nl),
    `25%` = rep(NA, nt * nl),
    `50%` = rep(NA, nt * nl),
    `75%` = rep(NA, nt * nl),
    `97.5%` = rep(NA, nt * nl),
    mu = rep(NA, nt * nl),
    sd = rep(NA, nt * nl),
    check.names = FALSE
  )
  r_estdf$time <- rep(seq_len(nt), each = nl)
  r_estdf$location <- rep(seq_len(nl), nt)
  for (time in seq_len(nt)) {
    for (location in seq_len(nl)) {
      r_estdf[r_estdf$time == time & r_estdf$location == location, 3:7] <- r_est[, time, location]
      r_estdf[r_estdf$time == time & r_estdf$location == location, "mu"] <- r_mu[time, location]
      r_estdf[r_estdf$time == time & r_estdf$location == location, "sd"] <- r_sd[time, location]
    }
  }
  r_estdf$param <- "R"
  epsilon_est <- quantile(
    fit$epsilon, probs = probs, na.rm = na.rm
  )
  eps_df <- data.frame(time = NA, location = NA)
  eps_df <- cbind(eps_df, epsilon_est)
  ## Stupid tall. make wide
  eps_df <- tibble::rownames_to_column(eps_df)
  eps_df <- tidyr::spread(eps_df, rowname, epsilon_est)
  eps_df$mu <- mean(fit$epsilon, na.rm  = na.rm)
  eps_df$sd <- sd(fit$epsilon, na.rm  = na.rm)
  eps_df$param <- "epsilon"
  rbind(eps_df, r_estdf)
}

summarise_vec <- function(vec, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) {
  vec_est <- quantile(vec, probs = probs, na.rm = na.rm)
  ## Tall. make wide
  eps_df <- tibble::rownames_to_column(data.frame(vec_est, check.names = FALSE))
  eps_df <- tidyr::spread(eps_df, rowname, vec_est)
  eps_df$mu <- mean(vec, na.rm  = na.rm)
  eps_df$sd <- sd(vec, na.rm  = na.rm)
  eps_df
}

summarise_epsilon <- function(fit, ...) {
  eps_df <- summarise_vec(fit$epsilon)
  eps_df$param <- "epsilon"
  eps_df
}

## Summarise epsilon - true_epsilon
## true_eps is the true epsilon value.
summarise_epsilon_error <- function(fit, true_eps, ...) {
  eps_df <- summarise_vec(fit$epsilon - true_eps)
  eps_df$param <- "epsilon_error"
  eps_df
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

  incid <- array(NA, dim = c(ndays, nlocations, nvariants))

  for (loc in seq_len(nlocations)) {
    for (v in seq_len(nvariants)) {
      incid[ ,loc, v] <- rbind(
        tail(incid_init[[v]]$counts, 1),
        as.matrix( #
          project(
            incid_init[[v]],
            ## R in the future so removing time of seeding
            R = rmatrix[-1, loc, v],
            si = simatrix[, v],
            n_sim = 1,
            n_days = ndays - 1,
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


#' Reorder an array of incidence data so that the most
#' transmissible variant is ordered first
#'
#' @param incidence An array containing time series of incidence
#' data with dimensions of 1) number of days; 2) number of locations;
#' 3) number of variants
#' @param t_start Start of time window over which to estimate
#' transmissibility (default is the minimum value = 2).
#' @param t_end End of time window over which to estimate transmissbility.
#' @si_distr serial interval distributions for the variants of interest. 1
#' column for each variant. si_distr[1,] must equal 0.
#'
#' @return array of dimensions ndays X nlocations X nvariants, where
#' the most transmissible variant in the time window of interest has
#' dimensions [,,1].

reorder_incidence <- function(incidence, t_start = 2, t_end, si_distr) {

  #identify variant with the largest estimated transmissibility over full time window
  R_init <- sapply(seq_len(dim(incidence)[3]), function(i) suppressWarnings(
    EpiEstim::estimate_R(apply(incidence[, , i, drop = FALSE], c(1, 3), sum)[,1],
                         method = "non_parametric_si",
                         config = EpiEstim::make_config(
                           si_distr = si_distr[, i],
                           t_start = t_start,
                           t_end = t_end)))$R$'Mean(R)')

  max_transmiss <- which.max(R_init)

  # reorder variants so most transmissible is first
  incid_reordered <- array(NA,
                           dim = dim(incidence))
  incid_reordered[, , 1] <- incidence[, , max_transmiss]
  incid_reordered[, , -1] <- incidence[, , -max_transmiss]

  incid_reordered

}
