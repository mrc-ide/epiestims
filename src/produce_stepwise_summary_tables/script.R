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
########### ONE LOCATION STEPWISE
#################################################
#################################################
#################################################
one_loc_step_eps <- readRDS("one_location_step_eps_summary_by_all_vars.rds")
one_loc_step_eps <- ungroup(one_loc_step_eps)
one_loc_step_eps <- filter(one_loc_step_eps, pt_est != "NaN")
one_loc_step_eps$rt_change <- paste(one_loc_step_eps$rt_ref,
                                    one_loc_step_eps$rt_post_step,
                                    sep = " -> ")

one_loc_step_eps <- select(one_loc_step_eps, rt_change, tmax, true_eps, pt_est:upper)
one_loc_step_eps <- mutate_if(one_loc_step_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
one_loc_step_pt <- select(one_loc_step_eps, rt_change, true_eps, tmax, pt_est)
one_loc_step_pt <- spread(one_loc_step_pt, tmax, pt_est)
one_loc_step_fill <- mutate_at(one_loc_step_pt, vars(`10`:`50`), f)

one_loc_step_eps$formatted <- pretty_ci(
  one_loc_step_eps$pt_est, one_loc_step_eps$lower, one_loc_step_eps$upper
)

x <- select(one_loc_step_eps, rt_change, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_change", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Rt change", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Rt change"]]))
y <- split(
  one_loc_step_fill, list(one_loc_step_fill$rt_change)
)

## Make dummy plot to extract legend
p <- ggplot(one_loc_step_eps) +
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

## Create tables
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo, col_start = 3)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/one_location_step_prop_in_95_1_legend.png",
          width = 540
        )
    }
    ggexport(
      tab,
      filename = glue("figures/one_location_step_prop_in_95_{i}.png"),
      width = 540
    )
  }
)


#################################################
#################################################
########### TWO LOCATION STEPWISE
#################################################
#################################################
#################################################
two_loc_step_eps <- readRDS("two_location_step_eps_summary_by_all_vars.rds")
two_loc_step_eps <- ungroup(two_loc_step_eps)
two_loc_step_eps <- filter(two_loc_step_eps, pt_est != "NaN")
two_loc_step_eps$rt_change <- paste(two_loc_step_eps$rt_ref_l1,
                                    two_loc_step_eps$rt_post_step_l1,
                                    sep = " -> ")

two_loc_step_eps <- select(two_loc_step_eps, rt_change, tmax, true_eps, pt_est:upper)
two_loc_step_eps <- mutate_if(two_loc_step_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
two_loc_step_pt <- select(two_loc_step_eps, rt_change, true_eps, tmax, pt_est)
two_loc_step_pt <- spread(two_loc_step_pt, tmax, pt_est)
two_loc_step_fill <- mutate_at(two_loc_step_pt, vars(`10`:`50`), f)

two_loc_step_eps$formatted <- pretty_ci(
  two_loc_step_eps$pt_est, two_loc_step_eps$lower, two_loc_step_eps$upper
)

x <- select(two_loc_step_eps, rt_change, true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_change", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Rt change", "True advantage",
                 "10", "20", "30", "40", "50")

x <- split(x, list(x[["Rt change"]]))
y <- split(
  two_loc_step_fill, list(two_loc_step_fill$rt_change)
)

## Create tables
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo, col_start = 3)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/two_location_step_prop_in_95_1_legend.png",
          width = 540
        )
    }
    ggexport(
      tab,
      filename = glue("figures/two_location_step_prop_in_95_{i}.png"),
      width = 540
    )
  }
)


#################################################
#################################################
########### TWO LOCATION STEPWISE (DIFFERENT)
#################################################
#################################################
#################################################
two_loc_step_diff_eps <- readRDS("two_location_step_diff_eps_summary_by_all_vars.rds")
two_loc_step_diff_eps <- ungroup(two_loc_step_diff_eps)
two_loc_step_diff_eps <- filter(two_loc_step_diff_eps, pt_est != "NaN")

two_loc_step_diff_eps$rt_change_loc1 <- paste(two_loc_step_diff_eps$rt_ref_l1,
                                              two_loc_step_diff_eps$rt_post_step_l1,
                                              sep = " -> ")
two_loc_step_diff_eps$rt_change_loc2 <- paste(two_loc_step_diff_eps$rt_ref_l2,
                                              two_loc_step_diff_eps$rt_post_step_l2,
                                              sep = " -> ")


two_loc_step_diff_eps <- select(two_loc_step_diff_eps, rt_change_loc1, rt_change_loc2,
                                tmax, true_eps, pt_est:upper)
two_loc_step_diff_eps <- mutate_if(two_loc_step_diff_eps, is.numeric, round, round_to)

## First spread, then pick color for each cell.
two_loc_step_diff_pt <- select(two_loc_step_diff_eps, rt_change_loc1, rt_change_loc2,
                               true_eps, tmax, pt_est)
two_loc_step_diff_pt <- spread(two_loc_step_diff_pt, tmax, pt_est)
two_loc_step_diff_fill <- mutate_at(two_loc_step_diff_pt, vars(`10`:`50`), f)

two_loc_step_diff_eps$formatted <- pretty_ci(
  two_loc_step_diff_eps$pt_est, two_loc_step_diff_eps$lower, two_loc_step_diff_eps$upper
)

x <- select(two_loc_step_diff_eps, rt_change_loc1, rt_change_loc2,
            true_eps, tmax, formatted)
x <- pivot_wider(
  x, id_cols = c("rt_change_loc1", "rt_change_loc2", "true_eps"),
  names_from = "tmax", values_from = "formatted"
)
## Nice column names
colnames(x) <- c("Rt (location 1)", "Rt (location 2)",
                 "True advantage", "10", "20", "30", "40", "50")

x <- split(x, list(x[["Rt (location 1)"]], x[["Rt (location 2)"]]))
y <- split(
  two_loc_step_diff_fill, list(two_loc_step_diff_fill$rt_change_loc1,
                               two_loc_step_diff_fill$rt_change_loc2)
)

## Create tables
pwalk(
  list(df = x, fillinfo = y, i = seq_along(x)),
  function(df, fillinfo, i) {
    tab <- prop_in_ci_table(df, fillinfo, col_start = 4)
    if (i == 1) {
      ggarrange(
        legend, tab, ncol = 1, nrow = 2, heights = c(0.1, 1)
      ) %>%
        ggexport(
          filename = "figures/two_location_step_diff_prop_in_95_1_legend.png",
          width = 540
        )
    }
    ggexport(
      tab,
      filename = glue("figures/two_location_step_diff_prop_in_95_{i}.png"),
      width = 540
    )
  }
)



if (! is.null(dev.list())) dev.off()
