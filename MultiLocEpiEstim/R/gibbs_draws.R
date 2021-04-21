default_priors <- function() {
  list(epsilon = list(shape = 1, scale = 1),
       Rt = list(shape = 1, scale = 1))
}

compute_lambda <- function(incid, w) {
  lambda <- array(NA, dim = dim(incid))
  for(l in seq_len(dim(incid)[2])) {
    for(v in seq_len(dim(incid)[3])) {
      lambda[, l, v] <- EpiEstim::overall_infectivity(incid[, l, v], w[, v])
    }
  }
  lambda
}

draw_epsilon <- function(R, incid, lambda, priors) {
  shape <- priors$epsilon$shape + sum(incid[, , 1])
  rate <- 1 / priors$epsilon$scale +
    sum(R * lambda[, , 1, drop = FALSE], na.rm = TRUE)
  scale <- 1 / rate
  rgamma(1, shape = shape, scale = scale)
}




n_v <- 2
n_loc <- 3
T <- 20
priors <- default_priors()
incid <- array(rpois(n_v * n_loc * T, 10), dim = c(T, n_loc, n_v))
w_v <- c(0, 0.2, 0.5, 0.3)
for(i in seq(2, n_v)) {
  w <- cbind(w, w_v)
}
lambda <- compute_lambda(incid, w)
R <- array(runif((n_v - 1) * n_loc * T, 0.5, 2.5), dim = c(T, n_loc, n_v - 1))
R[1, , ] <- NA # no estimates of R on first time step

draw_epsilon(R, incid, lambda, priors)
