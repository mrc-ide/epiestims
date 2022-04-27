## orderly::orderly_develop_start(use_draft = "newer")
source("mv_epiestim_params.R")
epi_params <- readRDS('Epi_param.rds')

infiles <- list(
  french = 'I_fr.rds',
  uk_alpha_wild = 'I_UK1.rds',
  uk_delta_alpha = 'I_UK2.rds'
)

incidence <- map(infiles, readRDS)

incid_array <- readRDS("incidence_array.rds")

## Analysis for each region
estimates <- map2(
  incid_array, incidence, function(x, df) {
    locations <- names(df[[1]])[-1]
    out <- imap(
      locations, function(location, index) {
        message("Location ", location)
        estimate_advantage(
          incid = x[, index, , drop = FALSE],
          si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = t_min
        )
      }
    )
    names(out) <- locations
    out
  }
)

saveRDS(estimates, "epsilon_estimates_per_region.rds")

eps_qntls <- map2(
  estimates,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(region, variants) {
    map_dfr(region, function(x) {
      summarise_epsilon(x, variants)
    }, .id = "region")
  })


saveRDS(eps_qntls, "epsilon_qntls_per_region.rds")

## Also do it once for all regions pooled
## for easy compilation
all_regions <- map(
  incid_array, function(x) {
    estimate_advantage(
      incid = x,
      si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
      mcmc_control = mcmc_controls,
      priors = priors,
      t_min = t_min
    )
  }
)
saveRDS(all_regions, "epsilon_estimates_whole_country.rds")

eps_qntls <- map2(
  all_regions,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(x, variants) {
    summarise_epsilon(x, variants)
  })


saveRDS(eps_qntls, "epsilon_qntls_whole_country.rds")

### Analysis for each quarter
periods <- readRDS("periods.rds")

eps_quarte <- map2(
  incid_array, periods, function(incid, period) {
    t_min <- period$intervals[-1] ## Remove the first element
    t_min <- head(t_min, -1)
    t_max <- period$intervals[-c(1, 2)]

    out <- map2(
      t_min, t_max, function(tmin, tmax) {
        message("tmin = ", tmin)
        message("tmax = ", tmax)
        message("nrow(incid) = ", nrow(incid))
        estimate_advantage(
          incid = incid,
          si_distr = cbind_rep(x = epi_params$SI, n = dim(incid)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = as.integer(tmin),
          t_max = as.integer(tmax)
        )
      }
    )
    names(out) <- paste("Quarter", 1:4)
    out
 }
)

names(eps_quarte) <- names(incid_array)

eps_qntls <- map2(
  eps_quarte,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(x, variants) {
    map_dfr(x, function(qx) {
      out <- apply(
        qx[["epsilon"]], 1, quantile,
        prob = c(0.025, 0.5, 0.975)
      )
      out <- data.frame(out)
      names(out) <- variants
      out <- tibble::rownames_to_column(out, "qntl")
      out <- tidyr::gather(out, variant, epsilon, -qntl)
      spread(out, qntl, epsilon)
    }, .id = "time_period")
  })

saveRDS(eps_qntls, "epsilon_qntls_time_periods.rds")


eps_over_time <- map2(
  incid_array, incidence, function(x, df) {
    locations <- names(df[[1]])[-1]
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    out <- imap(
      locations,
      function(location, index) {
        res <- map(t_max, function(tmax) {
          message("Location ", location, " tmax = ", tmax)
          estimate_advantage(
            incid = x[, index, , drop = FALSE],
            si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
            mcmc_control = mcmc_controls,
            priors = priors,
            t_min = as.integer(t_min),
            t_max = as.integer(tmax)
          )
        }
        )
        names(res) <- t_max
        res
      }
    )
    names(out) <- locations
    out
 }
)
