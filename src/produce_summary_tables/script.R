## orderly::orderly_develop_start()
## Aesthetics
## df is a grouped dataframe with column med which is the
## median error
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
ms_si <- c("X 0.5", "X 1.2", "X 1.5")
vary_si_eps <- readRDS("vary_si_eps_summary_by_all_vars.rds")
eps_vals <- unique(vary_si_eps$true_eps)
vary_si_eps$true_eps <- factor(
  vary_si_eps$true_eps,
  levels = eps_vals, ordered = TRUE
)
vary_si_eps$rt_ref <- factor(
  vary_si_eps$rt_ref
)

vary_si_eps$label <- multiplier_label(
  vary_si_eps$si_mu_variant, si_mu_ref
)


vary_si_eps <- vary_si_eps[vary_si_eps$label %in% ms_si, ]
vary_si_eps <- ungroup(vary_si_eps)
vary_si_eps <- select(vary_si_eps, rt_ref, label, tmax, true_eps, pt_est:upper)
vary_si_eps <- mutate_if(vary_si_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
vary_si_pt <- select(vary_si_eps, rt_ref, label, true_eps, tmax, pt_est)
vary_si_pt <- spread(vary_si_pt, tmax, pt_est)
## Make a table of colors
f <- scales::col_numeric("Greens", domain = c(0, 1))
vary_si_fill <- mutate_at(vary_si_pt, vars(`10`:`50`), f)

vary_si_eps$formatted <- glue(
  "{vary_si_eps$pt_est} \n ({vary_si_eps$lower}, {vary_si_eps$upper})"
)

x <- select(vary_si_eps, rt_ref, label, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "label", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)

x <- split(x, list(x$rt_ref, x$label))
y <- split(
  vary_si_fill, list(vary_si_fill$rt_ref, vary_si_fill$label)
)

pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- ggtexttable(
      df, rows = NULL, theme = ttheme(base_size = 6),
      cols = c("Reference Rt", "Variant SI Mean", "True epsilon",
               "10", "20", "30", "40", "50")
    )
    for (row in seq_len(nrow(df))) {
      for (col in seq(4, ncol(df))) {
        fill <- as.character(fillinfo[row, col])
        ## row + 1 is needed because header is in fact row 1
        tab <- table_cell_bg(
          tab, row = row + 1, column = col, fill = fill,
          alpha = 0.1
        )
      }
    }
    ggsave(glue("vary_si_prop_in_95_{i}.png"), tab)
  }
)


## Separate the case when si_mu_variant = si_ref_variant
same_si_mu <- vary_si_eps[vary_si_eps$label == "X 1", ]
same_si_mu1 <- same_si_mu[same_si_mu$tmax == ms_tmax, ]
same_si_mu2 <- same_si_mu[same_si_mu$tmax != ms_tmax, ]


## Ordered factors are being converted to
## integers when dumping them into a file.
vary_si_eps_tab$true_eps <- eps_vals
cat(
  stargazer(
    vary_si_eps_tab[, c("true_eps", "10", "50")], summary = FALSE
  ),
  file = "vary_si_epsor.tex"
)

## Main text figures
vary_si_eps1 <- vary_si_eps[vary_si_eps$tmax == ms_tmax, ]
## summarise_median_eps(vary_si_eps1)
# A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
## 1     -0.064     -0.003       0.006
## Top panel, mean for both variant and widltype
## the same.
p1a <- true_epsilon_vs_epsor(same_si_mu1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    ) +
  theme(legend.position = "none")
save_multiple(p1a, "figures/same_si_epsor")
## Middle panel, variant SI mean = X wildtype SI mean
p1b <- true_epsilon_vs_epsor(vary_si_eps1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )

save_multiple(p1b, "figures/vary_si_epsor")

## Supplementary figures
## Show epsor reducing by tmax
psi <- true_epsilon_vs_epsor(vary_si_eps, "Variant SI Mean") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(psi, "figures/vary_si_epsor_by_tmax")

## Classification
vary_si_classified <- readRDS("vary_si_classified.rds")
p <- classification_fig(vary_si_classified)
save_multiple(p, "figures/vary_si_classification")
######################################################################
######################################################################
################## VARY OFFSPRING ####################################
######################################################################
vary_offs_eps <- readRDS("vary_offs_eps_summary_by_all_vars.rds")
vary_offs_eps$true_eps <- factor(
  vary_offs_eps$true_eps, levels = eps_vals, ordered = TRUE
)
vary_offs_eps$rt_ref <- factor(vary_offs_eps$rt_ref)
vary_offs_eps$label <- round(vary_offs_eps$kappa, 1)
vary_offs_eps$label <- factor(vary_offs_eps$label)

vary_offs_ms <- vary_offs_eps[vary_offs_eps$tmax == ms_tmax, ]
## summarise_median_eps(vary_offs_ms)
# A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
## 1     -0.323     -0.007       0.001

## group_by(vary_offs_ms, kappa) %>% summarise_median_eps
## # A tibble: 3 × 4
##   kappa median_low median_med median_high
##   <dbl>      <dbl>      <dbl>       <dbl>
## 1   0.1     -0.376     -0.018       0.001
## 2   0.5     -0.176     -0.004       0.001
## 3   1       -0.163     -0.004       0.001

vary_offs_si <- vary_offs_eps[vary_offs_eps$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_epsor(vary_offs_ms, "Over-dispersion") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/vary_si_offs")

## Epsor over tmax
p1b <- true_epsilon_vs_epsor(vary_offs_si, "Over-dispersion") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1a, "figures/vary_si_offs_by_tmax")

vary_offs_eps_tab <- select(vary_offs_eps, -label) %>%
  group_by(true_eps, tmax) %>%
 summarise_median_eps() %>%
  format_median_eps()

## Ordered factors are being converted to
## integers when dumping them into a file.
vary_offs_eps_tab$true_eps <- eps_vals
cat(
  stargazer(
    vary_offs_eps_tab[, c("true_eps", "10", "50")],
    summary = FALSE, row.names = FALSE
  ),
  file = "vary_offs_epsor.tex"
)
vary_offs_classified <- readRDS("vary_offs_classified.rds")
p <- classification_fig(vary_offs_classified)
save_multiple(p, "figures/vary_offs_classification")
######################################################################
######################################################################
################## VARY CV ############################################
######################################################################
ms_cv <- c("X 0.5", "X 1.2", "X 1.5")
vary_cv_eps <- readRDS("vary_cv_eps_summary_by_all_vars.rds")
vary_cv_eps$true_eps <- factor(
  vary_cv_eps$true_eps, levels = eps_vals, ordered = TRUE
)
vary_cv_eps$rt_ref <- factor(vary_cv_eps$rt_ref)
vary_cv_eps$label <- multiplier_label(
  vary_cv_eps$si_cv_variant, si_std_ref / si_mu_ref
)
vary_cv_eps <- vary_cv_eps[vary_cv_eps$label %in% ms_cv, ]
vary_cv_eps$label <- factor(vary_cv_eps$label)

vary_cv_ms <- vary_cv_eps[vary_cv_eps$tmax == ms_tmax, ]
## summarise_median_eps(vary_cv_ms)
## # A tibble: 1 × 3
##   median_low median_med median_high
##        <dbl>      <dbl>       <dbl>
##      -0.101     -0.003       0.013

vary_cv_si <- vary_cv_eps[vary_cv_eps$tmax != ms_tmax, ]
p1a <- true_epsilon_vs_epsor(vary_cv_ms, "SI CV") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )
save_multiple(p1a, "figures/vary_si_cv")

## Epsor over tmax
p1b <- true_epsilon_vs_epsor(vary_cv_si, "SI CV") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(p1a, "figures/vary_si_cv_by_tmax")

vary_cv_eps_tab <- select(vary_cv_eps, -label) %>%
  group_by(true_eps, tmax) %>%
 summarise_median_eps() %>%
  format_median_eps()

## Ordered factors are being converted to
## integers when dumping them into a file.
vary_cv_eps_tab$true_eps <- eps_vals
cat(
  stargazer(
    vary_cv_eps_tab[, c("true_eps", "10", "50")],
    summary = FALSE, row.names = FALSE
  ),
  file = "vary_cv_epsor.tex"
)
vary_cv_classified <- readRDS("vary_cv_classified.rds")
p <- classification_fig(vary_cv_classified)
save_multiple(p, "figures/vary_cv_classification")
