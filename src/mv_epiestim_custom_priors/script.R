source("mv_epiestim_params.R")
epi_params <- readRDS('Epi_param.rds')

infiles <- list(
  french = 'I_fr.rds', uk_alpha_wild = 'I_UK1.rds',
  uk_delta_alpha = 'I_UK2.rds',
  french_betagamma = 'I_fr.rds'
)

incidence <- map(infiles, readRDS)

incid_array <- readRDS("incidence_array.rds")
incid_array[["french_betagamma"]] <- incid_array[["french"]][, , c(1, 3), drop = TRUE]
incid_array[["french"]] <- incid_array[["french"]][, , c(1, 2), drop = TRUE]

national_est <- readRDS("epsilon_estimates_whole_country.rds")
national_mu <- map(
  national_est, function(x) {
    mu <- apply(x[["epsilon"]], 1, mean)
    gamma_mucv2shapescale(
      mu = mu, cv = 0.1
    )
  }
)

custom_eps_priors <- list(
  french = list(
    shape = national_mu[["french"]]$shape,
    scale = national_mu[["french"]]$scale[1]
  ),
  uk_alpha_wild = national_mu[["uk_alpha_wild"]],
  uk_delta_alpha = national_mu[["uk_delta_alpha"]],
  french_betagamma = list(
    shape = national_mu[["french"]]$shape,
    scale = national_mu[["french"]]$scale[2]
  )
)

saveRDS(
  custom_eps_priors, "custom_epsilon_priors.rds"
)

estimates <- pmap(
  list(
    x = incid_array, y = incidence,
    z = custom_eps_priors),
  function(x, y, z) {
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    new_priors <- priors
    new_priors$epsilon$shape <- z$shape
    new_priors$epsilon$scale <- z$scale
    out <- map(
      t_max, function(tmax) {
        message("t_max = ", tmax)
        estimate_joint(
          incid = x,
          si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = new_priors,
          t_min = t_min,
          t_max = as.integer(tmax)
        )
    })
    names(out) <- y[[1]][["date"]][t_max]
    out
  }
)

saveRDS(estimates, "epsilon_estimates_over_time.rds")

eps_estimates <- map2(
  estimates,
  list(
    french = "alpha_vs_wild",
    uk_alpha_wild = "alpha_vs_wild",
    uk_delta_alpha = "delta_vs_alpha",
    french_betagmma = "beta-gamma_vs_wild"
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




cuml_incid <- readRDS("cuml_incid_all_variants.rds")
cols <- c("date", "wildtype", "betagamma", "proportion_wildtype", "proportion_betagamma")
cuml_incid[["french_betagamma"]] <- cuml_incid[["french"]][ ,cols]

## x-axis is now the proportion of the variant cases
joined <- map2(
  eps_estimates, cuml_incid, function(x, y) {
    x$date <- as.Date(x$date)
    left_join(x, y, by = "date")
  }
)

saveRDS(
  joined, "epsilon_estimates_with_variant_proportion.rds"
)
