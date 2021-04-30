process_fit <- function(fit, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  
  ## Get rid of the first row because of NAs
  r_est <- apply(
    fit$R[-1, , ], c(1, 2), quantile,
    probs = probs, na.rm = TRUE
  )
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
    check.names = FALSE
  )
  r_estdf$time <- rep(seq_len(nt), each = nl)
  r_estdf$location <- rep(seq_len(nl), nt)
  for (time in seq_len(nt)) {
    for (location in seq_len(nl)) {
      r_estdf[r_estdf$time == time & r_estdf$location == location, 3:7] <- r_est[, time, location]
    }
  }
  r_estdf$param <- "R"
  epsilon_est <- quantile(
    fit$epsilon, probs = probs, na.rm = TRUE
  )
  eps_df <- data.frame(time = NA, location = NA)
  eps_df <- cbind(eps_df, epsilon_est)
  ## Stupid tall. make wide
  eps_df <- tibble::rownames_to_column(eps_df)
  eps_df <- tidyr::spread(eps_df, rowname, epsilon_est)
  eps_df$param <- "epsilon"
  rbind(eps_df, r_estdf)
}