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

pretty_ci <- function(val, low, high, round_to = 2) {
  f <- function(x) {
    format(round(x, round_to), nsmall = 2)
  }
  glue("{f(val)} \n ({f(low)}, {f(high)})")
}

dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
round_to <- 2 ## Number of digits to round to
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
f <- scales::col_numeric("Greens", domain = c(0.5, 1))
vary_si_fill <- mutate_at(vary_si_pt, vars(`10`:`50`), f)

vary_si_eps$formatted <- pretty_ci(
  vary_si_eps$pt_est, vary_si_eps$lower, vary_si_eps$upper
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
      df, rows = NULL, theme = ttheme(base_size = 8),
      cols = c("Reference Rt", "Variant SI Mean", "True epsilon",
               "10", "20", "30", "40", "50")
    )
    for (row in seq_len(nrow(df))) {
      for (col in seq(4, ncol(df))) {
        fill <- as.character(fillinfo[row, col])
        ## row + 1 is needed because header is in fact row 1
        tab <- table_cell_bg(
          tab, row = row + 1, column = col, fill = fill,
          color = fill, alpha = 0.7
        )
      }
    }
    ggsave(
      glue("figures/vary_si_prop_in_95_{i}.png"),
      tab
    )
  }
)
