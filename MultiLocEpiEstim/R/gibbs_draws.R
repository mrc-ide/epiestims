default_priors <- function() {
  list(epsilon = list(shape = 1, scale = 1),
       R = list(shape = 1, scale = 1))
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


draw_epsilon <- function(R, incid, lambda, priors,
                         t_min = 2, t_max = nrow(incid),
                         seed = NULL) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- sum(incid[t, , 1]) + priors$epsilon$shape
  rate <- sum(R[t, , ] * lambda[t, , 1]) + 1 / priors$epsilon$scale
  scale <- 1 / rate
  rgamma(1, shape = shape, scale = scale)
}


draw_R <- function(epsilon, incid, lambda, priors,
                   t_min = 2, t_max = nrow(incid),
                   seed = NULL) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- apply(incid[t, , ], c(1, 2), sum) + priors$R$shape
  rate <- lambda[t, , 1] + epsilon * lambda[t, , 2] + 1 / priors$R$scale
  scale <- 1 / rate
  rgamma(length(t), shape = shape, scale = scale)
}
