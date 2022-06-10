## orderly::orderly_develop_start(use_draft = "newer")
source("mv_epiestim_params.R")
max_attempts <- 3
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
    newname <- glue("proportion_{colnames(out)[1]}")
    names(res)[names(res) == "prop_variant"] <- newname
  }
  res
}

epi_params <- readRDS('Epi_param.rds')

infiles <- list(
  french = 'I_fr.rds', uk_alpha_wild = 'I_UK1.rds',
  uk_delta_alpha = 'I_UK2.rds'
)

incidence <- map(infiles, readRDS)

## National incidence
fr_total_incid <- data.frame(
  date = incidence[["french"]][["wild"]][["date"]],
  wildtype = apply(incidence[["french"]][["wild"]][, -1], 1, sum),
  alpha = apply(incidence[["french"]][["alpha"]][, -1], 1, sum),
  betagamma = apply(incidence[["french"]][["beta/gamma"]][, -1], 1, sum)
)
## This is used to calculate proportion, so need
## to read in unadjusted numbers
eng_noadj <- readRDS("england_na_not_adjusted.rds")
uk1_total_incid <- eng_noadj[["alpha"]]
uk2_total_incid <- eng_noadj[["delta"]]
uk2_total_incid <- select(
  uk2_total_incid, date, alpha, delta, unknown
)

incid_array <- map(
  incidence, function(x) {
    time <- nrow(x[[1]])
    nlocation <- ncol(x[[1]]) - 1 ## -1 for date
    nvariant <- length(x)
    incid <- array(
      NA, dim = c(time, nlocation, nvariant)
    )
    for(i in 1:nvariant) {
      incid[,,i] <- as.matrix(x[[i]][, -1])
    }
  incid
})

saveRDS(incid_array, "incidence_array.rds")

## inftvty <- map2(
##   incid_array, incidence, function(x, df) {
##     nlocation <- dim(x)[2]
##     nvariant <- dim(x)[3]
##     si <- epi_params$SI
##     map_dfr(
##       1:nlocation, function(loc) {
##         map_dfr(1:nvariant, function(var) {
##           incid <- x[, loc, var, drop = TRUE]
##           data.frame(
##             inftvty = overall_infectivity(incid, si),
##             location = colnames(df[[var]])[1 + loc],
##             variant = names(df)[var],
##             date = df[[var]][["date"]]
##           )
##         })
##       })
##   }
## )

## x <- inftvty[[2]]
## x$date <- as.Date(x$date)
## p <- ggplot(x, aes(date, inftvty, col = variant)) +
##   geom_line() +
##   facet_wrap(~location, scales = "free_y") +
##   theme_minimal() +
##   theme(
##     legend.position = "top",
##     legend.title = element_blank(),
##     axis.title.x = element_blank()
##   ) +
##   ylab("Overall Infectivity")

window <- 6
nonovl_estimates <- map2(
  incid_array, incidence, function(x, df) {
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    out <- map(
      t_max, function(tmax) {
        message("t_max = ", tmax)
        out2 <- estimate_advantage(
          incid = x,
          si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = as.integer(tmax - window), # t_min,
          t_max = as.integer(tmax)
        )
        attempt <- 1
        while (! out2[["convergence"]]) {
          message("Attempt ", attempt)
          message("Not yet converged")
          mcmc_controls <- lapply(
              mcmc_controls, function(x) x * 2L
          )
          out2 <- estimate_advantage(
            incid = x,
            si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
            mcmc_control = mcmc_controls,
            priors = priors,
            t_min = as.integer(tmax - window), # t_min,
            t_max = as.integer(tmax)
          )
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            message("Aborting after 3 attempts")
            ## return whatever you've got.
            out2
          }
        }
      })
    names(out) <- df[[1]][["date"]][t_max]
    out
  }
)

saveRDS(
  nonovl_estimates,
  "nonoverlapping_epsilon_estimates.rds"
)

nonovl_eps_estimates <- map2(
  nonovl_estimates,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(estimate, variants) {
    map_dfr(estimate, function(x) {
      summarise_epsilon(x, variants)
    }, .id = "date")
  })
## Hacky, but will do for now.
nonovl_eps_estimates <- map(
  nonovl_eps_estimates, function(x) {
    x$date_min <- as.Date(x$date) - window
    x
})

saveRDS(nonovl_eps_estimates, "nonoverlapping_epsilon_qntls.rds")

nonovl_prop_variant <- map2(
  nonovl_estimates,
  list(
    french = fr_total_incid,
    uk_alpha_wild = uk1_total_incid,
    uk_delta_alpha = uk2_total_incid
  ), function(est, incid) {
    date_end <- as.Date(names(est))
    date_start <- date_end - window
    map2(
      date_start, date_end,
      function(start, end) {
        window_prop_variant(incid, start, end)
      }
    )
  }
)

saveRDS(nonovl_prop_variant, "nonoverlapping_prop_variant.rds")

estimates <- map2(
  incid_array, incidence, function(x, df) {
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    out <- map(
      t_max, function(tmax) {
        message("t_max = ", tmax)
        out2 <- estimate_advantage(
          incid = x,
          si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = as.integer(t_min),
          t_max = as.integer(tmax)
        )
        attempt <- 1
        while (! out2[["convergence"]]) {
          message("Attempt ", attempt)
          message("Not yet converged")
          mcmc_controls <- lapply(
              mcmc_controls, function(x) x * 2L
          )
          out2 <- estimate_advantage(
            incid = x,
            si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
            mcmc_control = mcmc_controls,
            priors = priors,
            t_min = as.integer(t_min),
            t_max = as.integer(tmax)
          )
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            message("Aborting after 3 attempts")
            ## return whatever you've got.
            out2
          }
        }
    })
    names(out) <- df[[1]][["date"]][t_max]
    out
  }
)

saveRDS(estimates, "epsilon_estimates_over_time.rds")

eps_estimates <- map2(
  estimates,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(estimate, variants) {
    map_dfr(estimate, function(x) {
      summarise_epsilon(x, variants)
    }, .id = "date")
  }
)

eps_estimates <- map2(eps_estimates, incidence, function(x, df) {
  x$date_min <- df[[1]][["date"]][t_min]
  x
}
)

saveRDS(eps_estimates, "epsilon_qntls_over_time.rds")

cuml_incid <- map(
  list(
    french = fr_total_incid,
    uk_alpha_wild = uk1_total_incid,
    uk_delta_alpha = uk2_total_incid
  ), function(x) {
    window_prop_variant(x, min(x$date), max(x$date))
  }
)

saveRDS(cuml_incid, "cuml_incid_all_variants.rds")

## x-axis is now the proportion of the variant cases
joined <- map2(
  eps_estimates, cuml_incid, function(x, y) {
    x$date <- as.Date(x$date)
    left_join(x, y, by = "date")
  }
)

saveRDS(joined, "epsilon_estimates_with_variant_proportion.rds")
