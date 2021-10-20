cbind_rep <- function(x, n) {
  matrix(x, nrow = length(x), ncol = n, byrow = FALSE)
}
## Setting this to 10 and incrementing in weeks makes sure we hit
## 31st December for Alpha.
t_min <- 10L
## priors <- default_priors()
priors <- list(
  epsilon = list(shape = 1, scale = 1),
  R = list(shape = 0.04, scale = 25)
)

mcmc_controls <- list(
  n_iter = 20000L, burnin = 5000L,
  thin = 10L
)
