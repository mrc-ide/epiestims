## orderly::orderly_develop_start()
## Aesthetics
## df is a grouped dataframe with column med which is the
## median error
## col_start is the first column index which
## has proportion in 95% CrI.
prop_in_ci_table <- function(df, fillinfo, col_start = 4) {
  tab <- ggtexttable(
    df, rows = NULL, theme = ttheme(base_size = 9)
  )
  for (row in seq_len(nrow(df))) {
    for (col in seq(col_start, ncol(df))) {
      fill <- as.character(fillinfo[row, col])
      ## row + 1 is needed because header is in fact row 1
      tab <- table_cell_bg(
        tab, row = row + 1, column = col, fill = fill,
        color = fill, alpha = 0.7
      )
    }
  }
  tab
}

## Make a table of colors
f <- scales::col_numeric("RdYlBu", domain = c(0, 1))

source("R/fig_utils.R")
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
ms_si <- c("X 0.5", "X 1.5", "X 2")
vary_si_eps <- readRDS("vary_si_eps_summary_by_all_vars.rds")
eps_vals <- unique(vary_si_eps$true_eps)
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

vary_si_fill <- mutate_at(vary_si_pt, vars(`10`:`50`), f)

vary_si_eps$formatted <- pretty_ci(
  vary_si_eps$pt_est, vary_si_eps$lower, vary_si_eps$upper
)

x <- select(vary_si_eps, rt_ref, label, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "label", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Variant SI Mean", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Variant SI Mean"]]))
y <- split(
  vary_si_fill, list(vary_si_fill$rt_ref, vary_si_fill$label)
)

## Make dummy plot to extract legend
p <- ggplot(vary_si_eps) +
  geom_tile(aes(tmax, true_eps, fill = pt_est)) +
  scale_fill_distiller(
    palette = "RdYlBu", limits = c(0, 1),
    direction = 1,
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1),
    name = "Proportion in 95% CrI"
  ) +
  theme_manuscript()

legend <- get_legend(p)

##ggsave("figures/legend.png", as_ggplot(legend))

pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/vary_si_prop_in_95_1_legend.png",
          width = 520
        )
    }
    ggexport(
      tab,
      filename = glue("figures/vary_si_prop_in_95_{i}.png"),
      width = 520
    )
  }
)

######################################################################
######################################################################
######################################################################
########### WRONG SI
######################################################################
wrong_si_eps <- readRDS("wrong_si_eps_summary_by_all_vars.rds")
eps_vals <- unique(wrong_si_eps$true_eps)
wrong_si_eps$label <- multiplier_label(
  wrong_si_eps$si_mu_variant, si_mu_ref
)


wrong_si_eps <- wrong_si_eps[wrong_si_eps$label %in% ms_si, ]
wrong_si_eps <- ungroup(wrong_si_eps)
wrong_si_eps <- select(wrong_si_eps, rt_ref, label, tmax, true_eps, pt_est:upper)
wrong_si_eps <- mutate_if(wrong_si_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
wrong_si_pt <- select(wrong_si_eps, rt_ref, label, true_eps, tmax, pt_est)
wrong_si_pt <- spread(wrong_si_pt, tmax, pt_est)

wrong_si_fill <- mutate_at(wrong_si_pt, vars(`10`:`50`), f)

wrong_si_eps$formatted <- pretty_ci(
  wrong_si_eps$pt_est, wrong_si_eps$lower, wrong_si_eps$upper
)

x <- select(wrong_si_eps, rt_ref, label, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "label", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Variant SI Mean", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Variant SI Mean"]]))
y <- split(
  wrong_si_fill, list(wrong_si_fill$rt_ref, wrong_si_fill$label)
)
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/wrong_si_prop_in_95_1_legend.png",
          width = 520
        )
    }
    ggexport(
      tab,
      filename = glue("figures/wrong_si_prop_in_95_{i}.png"),
      width = 520
    )
  }
)
#################################################
#################################################
########### VARY OFFSPRING
#################################################
#################################################
#################################################
vary_offs_eps <- readRDS("vary_offs_eps_summary_by_all_vars.rds")
vary_offs_eps <- ungroup(vary_offs_eps)
vary_offs_eps <- select(vary_offs_eps, rt_ref, kappa, tmax, true_eps, pt_est:upper)
vary_offs_eps <- mutate_if(vary_offs_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
vary_offs_pt <- select(vary_offs_eps, rt_ref, kappa, true_eps, tmax, pt_est)
vary_offs_pt <- spread(vary_offs_pt, tmax, pt_est)
vary_offs_fill <- mutate_at(vary_offs_pt, vars(`10`:`50`), f)

vary_offs_eps$formatted <- pretty_ci(
  vary_offs_eps$pt_est, vary_offs_eps$lower, vary_offs_eps$upper
)

x <- select(vary_offs_eps, rt_ref, kappa, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "kappa", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Overdispersion", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Overdispersion"]]))
y <- split(
  vary_offs_fill, list(vary_offs_fill$rt_ref, vary_offs_fill$kappa)
)
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/vary_offs_prop_in_95_1_legend.png",
          width = 520
        )
    }
    ggexport(
      tab,
      filename = glue("figures/vary_offs_prop_in_95_{i}.png"),
      width = 520
    )
  }
)

#################################################
#################################################
########### VARY CV
#################################################
#################################################
#################################################
ms_cv <- c("X 0.5", "X 1.5", "X 2")
vary_cv_eps <- readRDS("vary_cv_eps_summary_by_all_vars.rds")
vary_cv_eps <- ungroup(vary_cv_eps)
vary_cv_eps <- select(vary_cv_eps, rt_ref, si_cv_variant, tmax, true_eps, pt_est:upper)
vary_cv_eps$si_cv_variant <- multiplier_label(
  vary_cv_eps$si_cv_variant, si_std_ref / si_mu_ref
)
vary_cv_eps <- vary_cv_eps[vary_cv_eps$si_cv_variant %in% ms_cv, ]
vary_cv_eps <- mutate_if(vary_cv_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
vary_cv_pt <- select(vary_cv_eps, rt_ref, si_cv_variant, true_eps, tmax, pt_est)
vary_cv_pt <- spread(vary_cv_pt, tmax, pt_est)
vary_cv_fill <- mutate_at(vary_cv_pt, vars(`10`:`50`), f)

vary_cv_eps$formatted <- pretty_ci(
  vary_cv_eps$pt_est, vary_cv_eps$lower, vary_cv_eps$upper
)

x <- select(vary_cv_eps, rt_ref, si_cv_variant, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "si_cv_variant", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Variant CV", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Variant CV"]]))
y <- split(
  vary_cv_fill, list(vary_cv_fill$rt_ref, vary_cv_fill$si_cv_variant)
)
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/vary_cv_prop_in_95_1_legend.png",
          width = 520
        )
    }
    ggexport(
      tab,
      filename = glue("figures/vary_cv_prop_in_95_{i}.png"),
      width = 520
    )
  }
)
#################################################
#################################################
########### WRONG CV
#################################################
#################################################
#################################################
wrong_cv_eps <- readRDS("wrong_cv_eps_summary_by_all_vars.rds")
wrong_cv_eps <- ungroup(wrong_cv_eps)
wrong_cv_eps <- select(wrong_cv_eps, rt_ref, si_cv_variant, tmax, true_eps, pt_est:upper)
wrong_cv_eps$si_cv_variant <- multiplier_label(
  wrong_cv_eps$si_cv_variant, si_std_ref / si_mu_ref
)
wrong_cv_eps <- mutate_if(wrong_cv_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
wrong_cv_pt <- select(wrong_cv_eps, rt_ref, si_cv_variant, true_eps, tmax, pt_est)
wrong_cv_pt <- spread(wrong_cv_pt, tmax, pt_est)
wrong_cv_fill <- mutate_at(wrong_cv_pt, vars(`10`:`50`), f)

wrong_cv_eps$formatted <- pretty_ci(
  wrong_cv_eps$pt_est, wrong_cv_eps$lower, wrong_cv_eps$upper
)

x <- select(wrong_cv_eps, rt_ref, si_cv_variant, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "si_cv_variant", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Variant CV", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Variant CV"]]))
y <- split(
  wrong_cv_fill, list(wrong_cv_fill$rt_ref, wrong_cv_fill$si_cv_variant)
)
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/wrong_cv_prop_in_95_1_legend.png",
          width = 520
        )
    }
    ggexport(
      tab,
      filename = glue("figures/wrong_cv_prop_in_95_{i}.png"),
      width = 520
    )
  }
)
#################################################
#################################################
########### UNDERREPORTING
#################################################
#################################################
#################################################
underrep_eps <- readRDS("underrep_eps_summary_by_all_vars.rds")
underrep_eps <- ungroup(underrep_eps)
underrep_eps <- select(underrep_eps, rt_ref, p_report, tmax, true_eps, pt_est:upper)
underrep_eps$p_report <- round(
  underrep_eps$p_report, 1
)
underrep_eps <- mutate_if(underrep_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
underrep_pt <- select(underrep_eps, rt_ref, p_report, true_eps, tmax, pt_est)
underrep_pt <- spread(underrep_pt, tmax, pt_est)
underrep_fill <- mutate_at(underrep_pt, vars(`10`:`50`), f)

underrep_eps$formatted <- pretty_ci(
  underrep_eps$pt_est, underrep_eps$lower, underrep_eps$upper
)

x <- select(underrep_eps, rt_ref, p_report, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_ref", "p_report", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Reference Rt", "Reporting probability", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Reference Rt"]], x[["Reporting probability"]]))
y <- split(
  underrep_fill, list(underrep_fill$rt_ref, underrep_fill$p_report)
)
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/underrep_prop_in_95_1_legend.png",
          width = 540
        )
    }
    ggexport(
      tab,
      filename = glue("figures/underrep_prop_in_95_{i}.png"),
      width = 540
    )
  }
)

if (! is.null(dev.list())) dev.off()
