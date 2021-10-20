source("mv_epiestim_params.R")
epi_params <- readRDS('Epi_param.rds')

infiles <- list(
  french = 'I_fr.rds', uk_alpha_wild = 'I_UK1.rds',
  uk_delta_alpha = 'I_UK2.rds'
)

incidence <- map(infiles, readRDS)

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

estimates <- map2(
  incid_array, incidence, function(x, df) {
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    out <- map(
      t_max, function(tmax) {
        message("t_max = ", tmax)
        estimate_joint(
          incid = x,
          si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = t_min,
          t_max = as.integer(tmax)
        )
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
  function(region, variants) {
    map_dfr(region, function(x) {
      out <- apply(
        x[["epsilon"]], 1, quantile,
        prob = c(0.025, 0.5, 0.975)
      )
      out <- data.frame(out)
      names(out) <- variants
      out <- tibble::rownames_to_column(out, "qntl")
      out <- gather(out, variant, epsilon, -qntl)
      spread(out, qntl, epsilon)
    }, .id = "date")
  })

saveRDS(eps_estimates, "epsilon_qntls_over_time.rds")

## National incidence
fr_total_incid <- data.frame(
  date = incidence[["french"]][["wild"]][["date"]],
  wildtype = apply(incidence[["french"]][["wild"]][, -1], 1, sum),
  alpha = apply(incidence[["french"]][["alpha"]][, -1], 1, sum),
  betagamma = apply(incidence[["french"]][["beta/gamma"]][, -1], 1, sum)
)

uk1_total_incid <- data.frame(
  date = incidence[["uk_alpha_wild"]][["wild"]][["date"]],
  wildtype = apply(incidence[["uk_alpha_wild"]][["wild"]][, -1], 1, sum),
  alpha = apply(incidence[["uk_alpha_wild"]][["alpha"]][, -1], 1, sum)
)

uk2_total_incid <- data.frame(
  date = incidence[["uk_delta_alpha"]][["alpha"]][["date"]],
  alpha = apply(incidence[["uk_delta_alpha"]][["alpha"]][, -1], 1, sum),
  delta = apply(incidence[["uk_delta_alpha"]][["delta"]][, -1], 1, sum)
)



cuml_incid <- map(
  list(
    french = fr_total_incid,
    uk_alpha_wild = uk1_total_incid,
    uk_delta_alpha = uk2_total_incid
  ), function(x) {
    out <- x[, -1]
    ## Divide the daily incidence of variant
    ## by the daily incidence of wildtype
    prop_variant <- out / apply(out, 1, sum)
    res <- cbind(x, prop_variant)
    names(res) <- c(
      names(x), glue("proportion_{names(x[ ,-1])}")
    )
    res
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
