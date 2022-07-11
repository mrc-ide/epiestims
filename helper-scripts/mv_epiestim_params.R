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
## Proportion of variant in the time frame used for estimation
window_prop_variant <- function(incid, date_start, date_end) {
  x <- incid[incid$date >= date_start, ]
  x <- x[x$date <= date_end, ]
  out <- apply(x[, -1], 2, cumsum)
  res <- cbind(x, out)
  names(res) <- c(
    names(x), glue("cumulative_{names(x[, -1])}")
  )
  for (col in seq(2, ncol(out))) {
    ## Assume wildtype (or alpha, when estimating for delta)
    ## is always the first column.
    wt_plus_var <- out[, 1] + out[, col]
    res$prop_variant <- out[, col] / wt_plus_var
    newname <- glue("proportion_{colnames(out)[col]}")
    names(res)[names(res) == "prop_variant"] <- newname

    res$prop_variant <- out[, 1] / wt_plus_var
    newname <- glue("proportion_{colnames(out)[1]}_{colnames(out)[col]}")
    names(res)[names(res) == "prop_variant"] <- newname
  }
  res
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
  n_iter = 30000L, burnin = 7500L,
  thin = 10L
)
