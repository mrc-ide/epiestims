debug_compute_force_infection <- function (w, cases, R, t, instantaneous = FALSE)
{
  rev_w <- rev(w)
  ws <- utils::tail(rev_w, t)
  cases <- cases[seq_len(t), , drop = FALSE]
  R <- R[seq_len(t), , drop = FALSE]
  if(!instantaneous) {
    lambda <- ws %*% (cases * R)
  } else {
    lambda <- (ws %*% cases) * R[nrow(R)]
  }
  as.vector(lambda)
}

debug_project <- function (x, R, si, n_sim = 100, n_days = 7, R_fix_within = FALSE,
          model = c("poisson", "negbin"), size = 0.03, time_change = NULL,
          instantaneous_R = FALSE)
{
  model <- match.arg(model)
  if (!inherits(x, "incidence")) {
    msg <- "x is not an incidence object"
    stop(msg)
  }
  if (as.integer(mean(incidence::get_interval(x))) != 1L) {
    msg <- sprintf("daily incidence needed, but interval is %d days",
                   as.integer(mean(incidence::get_interval(x))))
    stop(msg)
  }
  if (ncol(incidence::get_counts(x)) > 1L) {
    msg <- sprintf("cannot use multiple groups in incidence object")
    stop(msg)
  }
  n_time_periods <- 1
  if (!is.null(time_change)) {
    if (!is.numeric(time_change)) {
      msg <- sprintf("`time_change` must be `numeric`, but is a `%s`",
                     paste(class(time_change), collapse = ", "))
      stop(msg)
    }
    n_time_periods <- length(time_change) + 1
    if (!is.vector(R)) {
      msg <- sprintf("`R` must be a `vector` or a `list` if `time_change` provided; it is a `%s`",
                     paste(class(R), collapse = ", "))
      stop(msg)
    }
    if (length(R) != n_time_periods) {
      msg <- sprintf("`R` must be a `list` of size %d to match %d time changes; found %d",
                     n_time_periods, n_time_periods - 1, length(R))
      stop(msg)
    }
  }
  projections:::assert_R(R)
  n_dates_x <- nrow(incidence::get_counts(x))
  t_max <- n_days + n_dates_x - 1
  if (inherits(si, "distcrete")) {
    if (as.integer(si$interval) != 1L) {
      msg <- sprintf("interval used in si is not 1 day, but %d)",
                     si$interval)
      stop(msg)
    }
    si <- si$d(1:t_max)
    si <- si/sum(si)
  } else {
    si <- si/sum(si)
    si <- c(si, rep(0, t_max - 1))
  }
  if (is.null(time_change)) {
    time_change <- Inf
  }
  I0 <- matrix(incidence::get_counts(x), nrow = n_dates_x,
               ncol = n_sim)
  out <- I0
  t_start <- n_dates_x + 1
  t_stop <- t_max + 1
  t_sim <- seq(from = t_start - 1, to = t_stop - 1, by = 1L)
  time_change <- t_start + time_change - 1
  if (!is.list(R)) {
    if (n_time_periods > 1L) {
      R <- as.list(R)
    } else {
      R <- list(R)
    }
  }
  #browser()
  if (all(is.finite(time_change))) {
    time_change_boundaries <- c(1, time_change, t_stop +
                                  1)
  } else {
    time_change_boundaries <- c(1, t_stop + 1)
  }
  #browser()
  R_t <- matrix(nrow = 0, ncol = n_sim)
  if (R_fix_within) {
    for (time_period in 1:n_time_periods) {
      R_time_period <- projections:::sample_(R[[time_period]], n_sim,
                               replace = TRUE)
      period_duration <- time_change_boundaries[time_period +
                                                  1] - time_change_boundaries[time_period]
      current_R_t <- do.call("rbind", replicate(period_duration,
                                                R_time_period, simplify = FALSE))
      R_t <- rbind(R_t, current_R_t)
    }
  } else {
    time_period <- 1L
    for (i in 1:t_stop) {
      R_time_period <- R[[time_period]]
      current_R_t <- projections:::sample_(R_time_period, n_sim, replace = TRUE)
      R_t <- rbind(R_t, current_R_t)
      if (i %in% time_change) {
        time_period <- time_period + 1
      }
    }
  }
  rownames(R_t) <- NULL
  for (i in t_sim) {
    lambda <- debug_compute_force_infection(si, out, R_t, i, instantaneous = instantaneous_R)
    # lambda2 <- EpiEstim::overall_infectivity(c(as.numeric(out[, 1]), 0), c(0, si)) * 3.75
    # lambda2 <- lambda2[length(lambda2)]
    # lambda3 <- EpiEstim::overall_infectivity(c(as.numeric(out[, 1]), 0), c(0, si)) * 2.25
    # lambda3 <- lambda3[length(lambda3)]
    if (model == "poisson") {
      out <- rbind(out, stats::rpois(n_sim, lambda))
    } else {
      size_adj <- lambda * size
      idx <- which(lambda == 0)
      size_adj[idx] <- 1
      out <- rbind(out, stats::rnbinom(n_sim, size = size_adj,
                                       mu = lambda))
    }
  }
  out <- out[(n_dates_x + 1):(n_dates_x + n_days), , drop = FALSE]
  dates <- utils::tail(incidence:::get_dates(x), 1) + seq_len(nrow(out))
  projections:::build_projections(out, dates, FALSE)
}


