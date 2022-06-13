cbind_rep <- function(x, n) {
  matrix(x, nrow = length(x), ncol = n, byrow = FALSE)
}

summarise_epsilon <- function(x, names) {
  out <- apply(
    x[["epsilon"]], 1, quantile,
    prob = c(0.025, 0.5, 0.975)
  )
  out <- data.frame(out)
  names(out) <- names
  out <- rownames_to_column(out, "qntl")
  out <- gather(out, variant, epsilon, -qntl)
  spread(out, qntl, epsilon)

}


## Setting this to 10 and incrementing in weeks makes sure we hit
## 31st December for Alpha.
t_min <- 10L
## priors <- default_priors()
priors <- list(
  epsilon = list(shape = 1, scale = 1),
  R = list(shape = 1, scale = 1)
)

mcmc_controls <- list(
  n_iter = 20000L, burnin = 5000L,
  thin = 10L
)
