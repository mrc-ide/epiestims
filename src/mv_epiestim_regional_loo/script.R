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

eps_non_overlapping <- map2(
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
          message("Leaving out Location ", location, " tmax = ", tmax)
          estimate_advantage(
            incid = x[, -index, , drop = FALSE],
            si_distr = cbind_rep(x = epi_params$SI, n = dim(x)[3]),
            mcmc_control = mcmc_controls,
            priors = priors,
            t_min = as.integer(tmax - 6),
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

saveRDS(eps_non_overlapping, "nonoverlapping_weekly_loo_epsilon.rds")

eps_non_overlapping_qntls <- map2(
  eps_non_overlapping,   list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
    ), function(x, variants) {
    map_dfr(x, function(y) {
      map_dfr(
        y, ~ summarise_epsilon(., variants),
        .id = "tmax"
      )
    }, .id = "leave_out")
  }
)

saveRDS(eps_non_overlapping_qntls, "nonoverlapping_weekly_loo_epsilon_qntls.rds")

## x <- eps_non_overlapping_qntls[[2]]
## x$tmax <- as.integer(x$tmax)

## all_regions <- readRDS("nonoverlapping_epsilon_qntls.rds")
## y <- all_regions[[2]]
## y$tmax <- seq(17, length.out = nrow(y), by = 7)
## y$leave_out <- "None"
## y <- y[, colnames(x)]
## x <- rbind(x, y)

## cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
##                 "#0072B2", "#D55E00", "#CC79A7")
## regions <- c("None", "East of England", "London", "Midlands", "North East and Yorkshire",
##   "North West", "South East", "South West")
## names(cbbPalette) <- regions

## x$leave_out <- factor(x$leave_out, levels = regions, ordered = TRUE)

## ggplot(x) +
##   geom_point(
##     aes(x = tmax, y = `50%`, col = leave_out),
##     position = position_dodge(width = 6)
##   ) +
##   geom_linerange(
##     aes(x = tmax, ymin = `2.5%`, ymax = `97.5%`, col = leave_out),
##     position = position_dodge(width = 6)
##   ) +
##   geom_hline(
##     yintercept = 1, linetype = "dashed", color = "red"
##   ) +
##   guides(
##     color = guide_legend(nrow = 1)
##   ) +
##   scale_color_manual(values = cbbPalette) +
##   ylab("Effective Transmission Advantage") +
##   xlab("tmax") +
##   theme_minimal() +
##   theme(legend.position = "top", legend.title = element_blank())




