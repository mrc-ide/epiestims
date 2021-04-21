context("Gibbs samplers")

test_that("draw_epsilon produces expected results", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  w <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, w)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_R produces expected results", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  w <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, w)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid, lambda, priors))
  x_mean <- Reduce("+", x) / length(x)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  expect_true(max(abs(x_mean[-c(1, 2, 3), ] - 1)) < 0.05)
})
