summarise_R <- function(fit, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) {
  ## Get rid of the first row because of NAs
  r_est <- apply(
    fit$R[-1, , ], c(1, 2), quantile,
    probs = probs, na.rm = na.rm
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
    fit$epsilon,
    probs = probs, na.rm = na.rm
  )
  eps_df <- data.frame(time = NA, location = NA)
  eps_df <- cbind(eps_df, epsilon_est)
  ## Stupid tall. make wide
  eps_df <- tibble::rownames_to_column(eps_df)
  eps_df <- tidyr::spread(eps_df, rowname, epsilon_est)
  eps_df$mu <- mean(fit$epsilon, na.rm = na.rm)
  eps_df$sd <- sd(fit$epsilon, na.rm = na.rm)
  eps_df$param <- "epsilon"
  rbind(eps_df, r_estdf)
}

summarise_vec <- function(vec, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) {
  vec_est <- quantile(vec, probs = probs, na.rm = na.rm)
  ## Tall. make wide
  eps_df <- rownames_to_column(data.frame(vec_est, check.names = FALSE))
  eps_df <- spread(eps_df, rowname, vec_est)
  eps_df$mu <- mean(vec, na.rm = na.rm)
  eps_df$sd <- sd(vec, na.rm = na.rm)
  eps_df
}
## This function assumes a single epsilon being returned.
## To fix this, we need to change summarise_vec to use apply
## Not fixing this in case it breaks something else. Do not use this
## function if you have more than one epsilon.
summarise_epsilon <- function(fit, ...) {
  eps_df <- summarise_vec(fit$epsilon)
  eps_df$param <- "epsilon"
  eps_df
}

## Summarise epsilon - true_epsilon
## true_eps is the true epsilon value.
## error can be absolute or relative;
summarise_epsilon_error <- function(fit, true_eps, error = "absolute", ...) {
  if (error != "absolute") {
    eps_df <- summarise_vec((fit$epsilon - true_eps) / true_eps)
    eps_df$param <- "epsilon_rel_error"
  } else {
    eps_df <- summarise_vec(fit$epsilon - true_eps)
    eps_df$param <- "epsilon_error"
  }
  eps_df
}

## df is grouped df, output of group_by
summarise_sims <- function(df) {
  summarise(
    df,
    n = sum(true_eps >= `2.5%` & true_eps <= `97.5%`),
    total = n(),
    pt_est = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 1],
    lower = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 2],
    upper = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 3],
    n50 = sum(true_eps >= `25%` & true_eps <= `75%`),
    ##total50 = n(),
    pt_est50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 1],
    lower50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 2],
    upper50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 3]
  )
}

estimate_uncertain <- function(df) {
  df$cv <- df$sd / df$mu
  ## Decided to not worry about this anymore
  df$confidence <- "High"
  df
}
