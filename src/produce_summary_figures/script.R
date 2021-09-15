## orderly::orderly_develop_start()
## Aesthetics
## df is a grouped dataframe with column med which is the
## median error
source("R/fig_utils.R")
summarise_median_err <- function(df, round_to = 3) {
  x <- summarise(
    df, median_low = quantile(med, 0.025),
    median_med = quantile(med, 0.5),
    median_high = quantile(med, 0.975)
  ) %>% ungroup()
  x <- mutate_if(x, is.numeric, round, round_to)
  x
}
## df is the output of summarise_median_err
format_median_err <- function(df) {
  df$formatted <-
  glue("{df$median_med}",
       " ({df$median_low}, {df$median_high})")
  df <- select(df, true_eps, tmax, formatted) %>%
    spread(key = tmax, value = formatted)
  df
}

dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
round_to <- 3 ## Number of digits to round to
#################################################
#################################################
######### SI MEAN SENSITIVITY ###################
#################################################
#################################################
## SIs of interest
ms_si <- c("X 0.5", "X 1.5", "X 2")
vary_si_err <- readRDS("vary_si_err_summary_by_all_vars.rds")
eps_vals <- unique(vary_si_err$true_eps)
vary_si_err$true_eps <- factor(
  vary_si_err$true_eps,
  levels = eps_vals, ordered = TRUE
)
vary_si_err$rt_ref <- factor(
  vary_si_err$rt_ref
)

vary_si_err$label <- multiplier_label(
  vary_si_err$si_mu_variant, si_mu_ref
)

## Separate the case when si_mu_variant = si_ref_variant
same_si_mu <- vary_si_err[vary_si_err$label == "X 1", ]
same_si_mu1 <- same_si_mu[same_si_mu$tmax == ms_tmax, ]
same_si_mu2 <- same_si_mu[same_si_mu$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_error(same_si_mu1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/same_si_error")

psi <- true_epsilon_vs_error(same_si_mu2, "Variant SI Mean") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(psi, "figures/same_si_error_by_tmax")

## At tmax = 50, across all epsilon values
## summarise_median_err(same_si_mu1)
## # A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
## 1     -0.054     -0.003       0.006
## These are the only ones we want to show, in Main and
## supplementary text
vary_si_err <- vary_si_err[vary_si_err$label %in% ms_si, ]
## Main text figures
vary_si_err1 <- vary_si_err[vary_si_err$tmax == ms_tmax, ]

## Variant SI mean = X wildtype SI mean
p1b <- true_epsilon_vs_error(vary_si_err1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )

save_multiple(p1b, "figures/vary_si_error")

## Supplementary figures
## Show error reducing by tmax
psi <- true_epsilon_vs_error(vary_si_err, "Variant SI Mean") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(psi, "figures/vary_si_error_by_tmax")

## Classification
vary_si_classified <- readRDS("vary_si_classified.rds")
##vary_si_classified <- vary_si_classified[vary_si_classified$confidence != "Low", ]
p <- classification_fig(vary_si_classified)
save_multiple(p, "figures/vary_si_classification")

######################################################################
######################################################################
################## WRONG SI ####################################
######################################################################
wrong_si_err <- readRDS("wrong_si_err_summary_by_all_vars.rds")
wrong_si_err$true_eps <- factor(
  wrong_si_err$true_eps, levels = eps_vals, ordered = TRUE
)
wrong_si_err$rt_ref <- factor(wrong_si_err$rt_ref)

wrong_si_err$label <- multiplier_label(
  wrong_si_err$si_mu_variant, si_mu_ref
)
wrong_si_err <- wrong_si_err[wrong_si_err$label %in% ms_si, ]
## Main text figures
wrong_si_err1 <- wrong_si_err[wrong_si_err$tmax == ms_tmax, ]
## Middle panel, variant SI mean = X wildtype SI mean
wrong_p1b <- true_epsilon_vs_error(wrong_si_err1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
      )
save_multiple(wrong_p1b, "figures/wrong_si_error")

## Common legend
wrong_p1b <- wrong_p1b +
  theme(axis.title.y = element_blank())
p <- ggarrange(p1b, wrong_p1b, ncol = 2, common.legend = TRUE)
save_multiple(plot = p, filename = "figures/correct_and_wrong_si_1")

## Faceted
wrong_si_err1$category <- "SI misspecified"
vary_si_err1$category <- "SI correctly specified"
vary_si_ms <- rbind(wrong_si_err1, vary_si_err1)
p <- true_epsilon_vs_error(vary_si_ms, "Variant SI Mean") +
      facet_grid(
      category~rt_ref, scales = "free_y",
      labeller = labeller(rt_ref = rt_labeller)
      )

save_multiple(plot = p, filename = "figures/correct_and_wrong_si")
## Supplementary figures
## Show error reducing by tmax
psi <- true_epsilon_vs_error(wrong_si_err, "Variant SI Mean") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(psi, "figures/wrong_si_error_by_tmax")

## Classification
wrong_si_classified <- readRDS("wrong_si_classified.rds")
p <- classification_fig(wrong_si_classified)
save_multiple(p, "figures/wrong_si_classification")
######################################################################
######################################################################
################## VARY OFFSPRING ####################################
######################################################################
vary_offs_err <- readRDS("vary_offs_err_summary_by_all_vars.rds")
vary_offs_err$true_eps <- factor(
  vary_offs_err$true_eps, levels = eps_vals, ordered = TRUE
)
vary_offs_err$rt_ref <- factor(vary_offs_err$rt_ref)
vary_offs_err$label <- round(vary_offs_err$kappa, 1)
vary_offs_err$label <- factor(vary_offs_err$label)

vary_offs_ms <- vary_offs_err[vary_offs_err$tmax == ms_tmax, ]
## summarise_median_err(vary_offs_ms)
# A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
## 1     -0.323     -0.007       0.001

## group_by(vary_offs_ms, kappa) %>% summarise_median_err
## # A tibble: 3 × 4
##   kappa median_low median_med median_high
##   <dbl>      <dbl>      <dbl>       <dbl>
## 1   0.1     -0.376     -0.018       0.001
## 2   0.5     -0.176     -0.004       0.001
## 3   1       -0.163     -0.004       0.001

vary_offs_si <- vary_offs_err[vary_offs_err$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_error(vary_offs_ms, "Over-dispersion") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/vary_si_offs")

## Error over tmax
p1b <- true_epsilon_vs_error(vary_offs_si, "Over-dispersion") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1b, "figures/vary_si_offs_by_tmax")

vary_offs_classified <- readRDS("vary_offs_classified.rds")
p <- classification_fig(vary_offs_classified)
save_multiple(p, "figures/vary_offs_classification")
######################################################################
######################################################################
################## VARY CV ############################################
######################################################################
ms_cv <- c("X 0.5", "X 1.5", "X 2")
vary_cv_err <- readRDS("vary_cv_err_summary_by_all_vars.rds")
vary_cv_err$true_eps <- factor(
  vary_cv_err$true_eps, levels = eps_vals, ordered = TRUE
)
vary_cv_err$rt_ref <- factor(vary_cv_err$rt_ref)
vary_cv_err$label <- multiplier_label(
  vary_cv_err$si_cv_variant, si_std_ref / si_mu_ref
)
vary_cv_err <- vary_cv_err[vary_cv_err$label %in% ms_cv, ]
vary_cv_err$label <- factor(vary_cv_err$label)
vary_cv_ms <- vary_cv_err[vary_cv_err$tmax == ms_tmax, ]
## summarise_median_err(vary_cv_ms)
## # A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
##      -0.101     -0.003       0.013

vary_cv_si <- vary_cv_err[vary_cv_err$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_error(vary_cv_ms, "SI CV") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/vary_si_cv")

## Error over tmax
p1b <- true_epsilon_vs_error(vary_cv_si, "SI CV") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1b, "figures/vary_si_cv_by_tmax")

vary_cv_classified <- readRDS("vary_cv_classified.rds")
vary_cv_classified$true_eps <- as.numeric(vary_cv_classified$true_eps)
p <- classification_fig(vary_cv_classified)
save_multiple(p, "figures/vary_cv_classification")
#################################################

wrong_cv_err <- readRDS("wrong_cv_err_summary_by_all_vars.rds")
wrong_cv_err$true_eps <- factor(
  wrong_cv_err$true_eps, levels = eps_vals, ordered = TRUE
)
wrong_cv_err$rt_ref <- factor(wrong_cv_err$rt_ref)
wrong_cv_err$label <- multiplier_label(
  wrong_cv_err$si_cv_variant, si_std_ref / si_mu_ref
)
wrong_cv_err <- wrong_cv_err[wrong_cv_err$label %in% ms_cv, ]
wrong_cv_err$label <- factor(wrong_cv_err$label)

wrong_cv_ms <- wrong_cv_err[wrong_cv_err$tmax == ms_tmax, ]
## summarise_median_err(wrong_cv_ms)
## # A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
##      -0.101     -0.003       0.013

wrong_cv_si <- wrong_cv_err[wrong_cv_err$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_error(wrong_cv_ms, "Variant SI CV") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/wrong_cv")

## Error over tmax
p1b <- true_epsilon_vs_error(wrong_cv_si, "Variant SI CV") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1b, "figures/wrong_cv_by_tmax")


wrong_cv_ms$category <- "SI CV misspecified"
vary_cv_ms$category <- "SI CV correctly specified"
vary_cv_ms <- rbind(wrong_cv_ms, vary_cv_ms)
p <- true_epsilon_vs_error(vary_cv_ms, "Variant SI CV") +
      facet_grid(
      category~rt_ref, scales = "free_y",
      labeller = labeller(rt_ref = rt_labeller)
      )
save_multiple(plot = p, filename = "figures/correct_and_wrong_cv")

wrong_cv_classified <- readRDS("wrong_cv_classified.rds")
wrong_cv_classified$true_eps <- as.numeric(wrong_cv_classified$true_eps)
p <- classification_fig(wrong_cv_classified)
save_multiple(p, "figures/wrong_cv_classification")
######################################################################
######################################################################
################## UNDERREPORTING ####################################
######################################################################
underrep_err <- readRDS("underrep_err_summary_by_all_vars.rds")
underrep_err$true_eps <- factor(
  underrep_err$true_eps, levels = eps_vals, ordered = TRUE
)
underrep_err$rt_ref <- factor(underrep_err$rt_ref)
underrep_err$label <- round(underrep_err$p_report, 1)
underrep_err$label <- factor(underrep_err$label)

## Error over tmax
p1b <- true_epsilon_vs_error(underrep_err, "Reporting probability") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1b, "figures/underrep_by_tmax")

underrep_classified <- readRDS("underrep_classified.rds")
p <- classification_fig(underrep_classified)
save_multiple(p, "figures/underrep_classification")

if (! is.null(dev.list())) dev.off()
