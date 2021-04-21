#' Set default for Gamma priors
#'
#' @return a list of default parameters for the priors.
#'   Values can then be manually be edited as in the examples below.
#'   Users could use functions `epitrix::gamma_shapescale2mucv` and
#'   `epitrix::gamma_mucv2shapescale` to set the shape and scale corresponding
#'   to the desired prior mean and coefficient of variation.
#' @export
#'
#' @examples
#' priors <- default_priors
#' # change the prior for R to have a mean of 3
#' priors$R$shape <- 3
#'
default_priors <- function() {
  list(epsilon = list(shape = 1, scale = 1),
       R = list(shape = 1, scale = 1))
}


#' Compute the overall infectivity
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param w a matrix with two columns, each containing the probability mass
#'   function for the discrete serial interval for each of the two
#'   pathogen/strain/variants, starting with the probability mass function
#'   for day 0 in the first row, which should be 0. each column in the matrix
#'   should sum to 1
#'
#' @return a multidimensional array containing values of the overall
#'   infectiousness for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectiousness for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectiousness of those
#'   past incident cases.
#'
#' @export
#'
#' @examples
compute_lambda <- function(incid, w) {
  ## TODO: check that w[1, ] == 0
  ## TODO: check that all(colSums(w) == 1)
  ## TODO: check that all(w >= 0)
  lambda <- array(NA, dim = dim(incid))
  for(l in seq_len(dim(incid)[2])) {
    for(v in seq_len(dim(incid)[3])) {
      lambda[, l, v] <- EpiEstim::overall_infectivity(incid[, l, v], w[, v])
    }
  }
  lambda
}


#' Draw epsilon from marginal posterior distribution
#'
#' @param R a matrix with dimensions containing values of the instantaneous
#'   reproduction number for each time step (row) and location (column), for
#'   the reference pathogen/strain/variant
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param lambda a multidimensional array containing values of the overall
#'   infectiousness for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectiousness for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectiousness of those
#'   past incident cases. It can be calculated from the incidence `incid` and
#'   the distribution of the serial interval using function `compute_lambda`)
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param t_min an integer >1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return a value of epsilon drawn from the marginal posterior distribution
#'
#' @export
#'
#' @examples
#'
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' w <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, w)
#' # Constant reproduction number of 1
#' R <- matrix(1, nrow = T, ncol = n_loc)
#' R[1, ] <- NA # no estimates of R on first time step
#' draw_epsilon(R, incid, lambda, priors, seed = 1)
#'
draw_epsilon <- function(R, incid, lambda, priors,
                         t_min = 2, t_max = nrow(incid),
                         seed = NULL) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  ## TODO: check R >=0
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- sum(incid[t, , 1]) + priors$epsilon$shape
  rate <- sum(R[t, ] * lambda[t, , 1]) + 1 / priors$epsilon$scale
  scale <- 1 / rate
  rgamma(1, shape = shape, scale = scale)
}


#' Draw R from marginal posterior distribution
#'
#' @param epsilon the value of the relative transmissibility of the "new"
#'   pathogen/strain/variant compared to the reference pathogen/strain/variant
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param lambda a multidimensional array containing values of the overall
#'   infectiousness for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectiousness for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectiousness of those
#'   past incident cases. It can be calculated from the incidence `incid` and
#'   the distribution of the serial interval using function `compute_lambda`)
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param t_min an integer >1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return a matrix of the instantaneous reproduction number R for the reference
#'   pathogen/strain/variant for each time step (row) and each location (column)
#'   drawn from the marginal posterior distribution
#'
#' @export
#'
#' @examples
#'
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' w <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, w)
#' # Epsilon = 1 i.e. no transmission advantage
#' epsilon <- 1
#' draw_R(epsilon, incid, lambda, priors, seed = 1)
#'
draw_R <- function(epsilon, incid, lambda, priors,
                   t_min = 2, t_max = nrow(incid),
                   seed = NULL) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  ## TODO: check epsilon >0
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- apply(incid[t, , ], c(1, 2), sum) + priors$R$shape
  shape_flat <- as.numeric(shape)
  rate <- lambda[t, , 1] + epsilon * lambda[t, , 2] + 1 / priors$R$scale
  scale <- 1 / rate
  scale_flat <- as.numeric(scale)
  R_flat <- rgamma(length(shape_flat), shape = shape_flat, scale = scale_flat)
  R_fill <- matrix(R_flat, nrow = nrow(shape), ncol = ncol(shape))
  R <- matrix(NA, nrow(incid), ncol(incid))
  R[t, ] <- R_fill
  R
}
