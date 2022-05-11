### Summary tables
naive_eps <- map(
  c(french = "naive_epsilon_fr.rds",
    uk_alpha_wild = "naive_epsilon_UK1.rds",
    uk_delta_alpha = "naive_epsilon_UK2.rds"), function(x) {
      readRDS(x)$sum_epsilon
    }
)

naive_eps[["french_betagamma"]] <- naive_eps[["french"]][["wild_beta/gamma"]]
naive_eps[["french"]] <- naive_eps[["french"]][["wild_alpha"]]
naive_eps[["uk_alpha_wild"]] <- naive_eps[["uk_alpha_wild"]][[1]]
naive_eps[["uk_delta_alpha"]] <- naive_eps[["uk_delta_alpha"]][[1]]

## MV-EpiEstim time periods Q1-4
mvepi_q <- readRDS("epsilon_qntls_time_periods.rds")
mvepi_q[["french_betagamma"]] <-
  mvepi_q[["french"]][mvepi_q[["french"]]$variant != "alpha_vs_wild", ]
mvepi_q[["french"]] <-
  mvepi_q[["french"]][mvepi_q[["french"]]$variant == "alpha_vs_wild", ]

pretty_ci <- function(val, low, high, round_to = 2) {
  f <- function(x) {
    format(round(x, round_to), nsmall = 2)
  }
  glue("{f(val)} ({f(low)}, {f(high)})")
}

## get dates corresponding to quarters
periods <- readRDS("periods.rds")
periods[["french_betagamma"]] <- periods[["periods_fr"]]

quarter_dates <- map2(
  incidence, periods, function(incid, period) {
    t_min <- period$intervals[-1] ## Remove the first element
    t_min <- head(t_min, -1)
    t_max <- period$intervals[-c(1, 2)]
    t_min <- format(incid$date[t_min], "%d-%b-%Y")
    t_max <- format(incid$date[t_max], "%d-%b-%Y")
    glue("{t_min} to {t_max}")
  }
)


si_tables <- pmap(
  list(naive = naive_eps, mvepi = regional, overall = national, qx = mvepi_q,
       qd = quarter_dates),
  function(naive, mvepi, overall, qx, qd) {
    overall$region <- "All"
    mvepi <- rbind(overall[, colnames(mvepi)], mvepi)
    qx$region <- qx$time_period
    mvepi <- rbind(qx[, colnames(mvepi)], mvepi)
    idx <- which(is.na(naive$med))
    naive$formatted <- pretty_ci(naive$med, naive$low, naive$up)
    naive$formatted[idx] <- "-"
    mvepi$formatted <- pretty_ci(mvepi$`50%`, mvepi$`2.5%`, mvepi$`97.5%`)
    naive <- select(naive, "Region/Time-Period" = name, `Naive` = formatted)
    mvepi <- select(mvepi, "Region/Time-Period" = region, `MV-EpiEstim` = formatted)
    out <- left_join(naive, mvepi, by = "Region/Time-Period")
    ## out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 1"] <- qd[1]
    ## out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 2"] <- qd[2]
    ## out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 3"] <- qd[3]
    ## out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 4"] <- qd[4]
    out
  }
)
