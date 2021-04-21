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
  # ignoring the first time step as we are conditioning on this
  shape <- sum(incid[-1, , 1]) + priors$epsilon$shape
  rate <- sum(R[-1, , ] * lambda[-1, , 1]) + 1 / priors$epsilon$scale
  scale <- 1 / rate
  rgamma(1, shape = shape, scale = scale)
}




n_v <- 2 # 2 variants
n_loc <- 3 # 3 locations
T <- 100 # 20 time steps

priors <- default_priors()

# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))

# arbitrary serial interval
w_v <- c(0, 0.2, 0.5, 0.3)
for(i in seq(2, n_v)) {
  w <- cbind(w, w_v)
}
lambda <- compute_lambda(incid, w)

# Reproduction number of 1
R <- array(1, dim = c(T, n_loc, n_v - 1))
R[1, , ] <- NA # no estimates of R on first time step

draw_epsilon(R, incid, lambda, priors)
x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))
mean(x)
median(x) ### This should be 1 - it is not 1 because of the first few timesteps AND because of priors - but close enough
