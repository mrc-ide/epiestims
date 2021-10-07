library(EpiEstim)
library(glue)
library(purrr)

t_min <- 14L
priors <- default_priors()
mcmc_controls <- list(
  n_iter = 20000L, burnin = 5000L, thin = 10L
)
epi_params <- readRDS('Rdata/Epi_param.rds')

infiles <- list(
  french = 'Rdata/I_fr.rds',
  uk_alpha_wild = 'Rdata/I_UK1.rds',
  uk_delta_alpha = 'Rdata/I_UK2.rds'
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
        estimate_joint(
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
      out <- apply(
        x[["epsilon"]], 1, quantile,
        prob = c(0.025, 0.5, 0.975)
      )
      out <- data.frame(out)
      names(out) <- variants
      out <- tibble::rownames_to_column(out, "qntl")
      out <- tidyr::gather(out, variant, epsilon, -qntl)
      tidyr::spread(out, qntl, epsilon)
    }, .id = "region")
  })


saveRDS(eps_qntls, "epsilon_qntls_per_region.rds")

