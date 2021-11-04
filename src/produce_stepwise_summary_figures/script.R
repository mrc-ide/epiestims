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

######################################################################
######################################################################
################## ONE LOCATION STEPWISE #############################
######################################################################
one_loc_step_err <- readRDS("one_loc_step_err_summary_by_all_vars.rds")
eps_vals <- unique(one_loc_step_err$true_eps)

one_loc_step_err$true_eps <- factor(
  one_loc_step_err$true_eps, levels = eps_vals, ordered = TRUE
)
one_loc_step_err$rt_change <- factor(paste(one_loc_step_err$rt_ref,
                                           one_loc_step_err$rt_post_step,
                                           sep = " to ")
)
one_loc_step_err$si_mu_variant <- 5.4
one_loc_step_err$label <- multiplier_label(
  one_loc_step_err$si_mu_variant, si_mu_ref
)
one_loc_step_err$metric <- "Bias" # variable for faceting plots

## Error over tmax and rt_change

psi <- true_epsilon_vs_error(one_loc_step_err, "Variant SI Mean") +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )
save_multiple(psi, "figures/one_loc_step_error_by_tmax")


## sd over tmax and rt change

one_loc_step_sd <- readRDS("one_loc_step_err_sd_summary_by_all_vars.rds")
eps_vals <- unique(one_loc_step_sd$true_eps)

one_loc_step_sd$true_eps <- factor(
  one_loc_step_sd$true_eps, levels = eps_vals, ordered = TRUE
)
one_loc_step_sd$rt_change <- factor(paste(one_loc_step_sd$rt_ref,
                                          one_loc_step_sd$rt_post_step,
                                           sep = " to ")
)
one_loc_step_sd$si_mu_variant <- 5.4
one_loc_step_sd$label <- multiplier_label(
  one_loc_step_sd$si_mu_variant, si_mu_ref
)
one_loc_step_sd$metric <- "Uncertainty"

psd <- true_epsilon_vs_sd(one_loc_step_sd) +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )
save_multiple(psd, "figures/one_loc_step_sd_by_tmax")

## Classification figures
one_loc_step_classified <- readRDS("one_loc_step_classified.rds")
one_loc_step_classified$rt_change <- factor(paste(one_loc_step_classified$rt_ref,
                                                  one_loc_step_classified$rt_post_step,
                                                  sep = " to ")
)
one_loc_step_classified$true_eps <- factor(
  one_loc_step_classified$true_eps, levels = eps_vals, ordered = TRUE
)
idx1 <- which(one_loc_step_classified$true_label == "No transmission advantage" &
                one_loc_step_classified$est_class == "Unclear")
idx2 <- which(one_loc_step_classified$true_label == one_loc_step_classified$est_class)
x <- one_loc_step_classified[c(idx1, idx2), ]
x$metric <- "Classification"
y <- split(x, x$rt_change)

iwalk(y, function(change, index) {
  
  p <- classification_fig(change)
  save_multiple(
    p, glue("figures/one_loc_step_classification_{index}")
  )
  
})


one_loc_step_classified <- rename(x,
                                med = PointEst, low = Lower, high = Upper)

# Coverage probability
one_loc_step_coverage <- readRDS("one_loc_step_eps_summary_by_all_vars.rds")
one_loc_step_coverage <- filter(one_loc_step_coverage, pt_est != "NaN")
one_loc_step_coverage$rt_change <- factor(paste(one_loc_step_coverage$rt_ref,
                                                  one_loc_step_coverage$rt_post_step,
                                                  sep = " to ")
)
one_loc_step_coverage$true_eps <- factor(
  one_loc_step_coverage$true_eps, levels = eps_vals, ordered = TRUE
)
one_loc_step_coverage$metric <- "Coverage probability"
one_loc_step_coverage <- rename(one_loc_step_coverage,
                                med = pt_est, low = lower, high = upper)


one_loc_step_coverage$rt_ref <- as.numeric(one_loc_step_coverage$rt_ref)
one_loc_step_coverage$rt_post_step <- as.numeric(one_loc_step_coverage$rt_post_step)


together <- rbind(one_loc_step_err, one_loc_step_sd,
                  one_loc_step_coverage, one_loc_step_classified)


together$metric <- factor(
  together$metric, levels = c("Bias", "Uncertainty",
                              "Coverage probability",
                              "Classification"),
  ordered = TRUE
)


dummy <- data.frame(
  metric = c("Bias", "Coverage probability"),
  ##true_eps = levels(y$true_eps),
  low = 0,
  high = 1
)
dummy2 <- data.frame(
  metric = c("Bias", "Coverage probability"),
  y = c(0, 0.95)
)
dummy$metric <- factor(
  dummy$metric, levels = levels(together$metric)
)
dummy2$metric <- factor(
  dummy2$metric, levels = levels(together$metric)
)



y <- split(together, list(together$tmax, together$rt_change))
iwalk(y, function(z, index) {

p <- ggplot(z) +
  geom_point(
    aes(true_eps, med), col = "black",
    position = position_dodge(width = dodge_width),
    size = 1.4
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high),
    col = "black",
    position = position_dodge(width = dodge_width),
    size = 1.1
  ) +
  geom_blank(
    data = dummy, aes(y = low)
  ) +
  geom_blank(
    data = dummy, aes(y = high)
  ) +
  geom_hline(
    data = dummy2, aes(yintercept = y),
    linetype = "dashed"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  theme_manuscript() +
  xlab("Transmission Advantage") +
  ylab("")

save_multiple(
  p, glue("figures/one_loc_step_panel_{index}")
)

})

######################################################################
######################################################################
################## TWO LOCATION STEPWISE #############################
######################################################################
two_loc_step_err <- readRDS("two_loc_step_err_summary_by_all_vars.rds")

two_loc_step_err$true_eps <- factor(
  two_loc_step_err$true_eps, levels = eps_vals, ordered = TRUE
)
two_loc_step_err$rt_change <- factor(paste(two_loc_step_err$rt_ref_l1,
                                           two_loc_step_err$rt_post_step_l1,
                                           sep = " to ")
)

two_loc_step_err$label <- multiplier_label(
  two_loc_step_err$si_mu_variant, si_mu_ref
)

## Error over tmax and rt_change

psi <- true_epsilon_vs_error(two_loc_step_err, "Variant SI Mean") +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )
save_multiple(psi, "figures/two_loc_step_error_by_tmax")


## sd over tmax and rt_change

two_loc_step_sd <- readRDS("two_loc_step_err_sd_summary_by_all_vars.rds")

two_loc_step_sd$true_eps <- factor(
  two_loc_step_sd$true_eps, levels = eps_vals, ordered = TRUE
)
two_loc_step_sd$rt_change <- factor(paste(two_loc_step_sd$rt_ref_l1,
                                          two_loc_step_sd$rt_post_step_l1,
                                           sep = " to ")
)

two_loc_step_sd$label <- multiplier_label(
  two_loc_step_sd$si_mu_variant, si_mu_ref
)

psd <- true_epsilon_vs_sd(two_loc_step_sd) +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )

save_multiple(psd, "figures/two_loc_step_sd_by_tmax")


## classification figure

two_loc_step_classified <- readRDS("two_loc_step_classified.rds")
two_loc_step_classified$rt_change <- factor(paste(two_loc_step_classified$rt_ref_l1,
                                                  two_loc_step_classified$rt_post_step_l1,
                                                  sep = " to ")
)
two_loc_step_classified$true_eps <- factor(
  two_loc_step_classified$true_eps, levels = eps_vals, ordered = TRUE
)
idx1 <- which(two_loc_step_classified$true_label == "No transmission advantage" &
                two_loc_step_classified$est_class == "Unclear")
idx2 <- which(two_loc_step_classified$true_label == two_loc_step_classified$est_class)
x <- two_loc_step_classified[c(idx1, idx2), ]
y <- split(x, x$rt_change)

iwalk(y, function(change, index) {
  
  p <- classification_fig(change)
  save_multiple(
    p, glue("figures/two_loc_step_classification_{index}")
  )
  
})

######################################################################
######################################################################
################## TWO LOCATION STEPWISE (DIFFERING LOCATIONS) #######
######################################################################
two_loc_step_diff_err <- readRDS("two_loc_step_diff_err_summary_by_all_vars.rds")

two_loc_step_diff_err$true_eps <- factor(
  two_loc_step_diff_err$true_eps, levels = eps_vals, ordered = TRUE
)

loc1_change <- paste(two_loc_step_diff_err$rt_ref_l1,
                     two_loc_step_diff_err$rt_post_step_l1,
                     sep = " -> ")
loc2_change <- paste(two_loc_step_diff_err$rt_ref_l2,
                     two_loc_step_diff_err$rt_post_step_l2,
                     sep = " -> ")


two_loc_step_diff_err$rt_change <- paste("Location 1: ", loc1_change,
                              ", Location 2: ", loc2_change)

two_loc_step_diff_err$label <- multiplier_label(
  two_loc_step_diff_err$si_mu_variant, si_mu_ref
)

## Error over tmax and rt_change

psi <- true_epsilon_vs_error(two_loc_step_diff_err, "Variant SI Mean") +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )
save_multiple(psi, "figures/two_loc_step_diff_error_by_tmax")


## sd over tmax and rt change

two_loc_step_diff_sd <- readRDS("two_loc_step_diff_err_sd_summary_by_all_vars.rds")

two_loc_step_diff_sd$true_eps <- factor(
  two_loc_step_diff_sd$true_eps, levels = eps_vals, ordered = TRUE
)

loc1_change <- paste(two_loc_step_diff_sd$rt_ref_l1,
                     two_loc_step_diff_sd$rt_post_step_l1,
                     sep = " -> ")
loc2_change <- paste(two_loc_step_diff_sd$rt_ref_l2,
                     two_loc_step_diff_sd$rt_post_step_l2,
                     sep = " -> ")


two_loc_step_diff_sd$rt_change <- paste("Location 1: ", loc1_change,
                                         ", Location 2: ", loc2_change)

two_loc_step_diff_sd$label <- multiplier_label(
  two_loc_step_diff_sd$si_mu_variant, si_mu_ref
)


psd <- true_epsilon_vs_sd(two_loc_step_diff_sd) +
  facet_grid(
    tmax~rt_change,
    labeller = labeller(rt_change = rt_change_labeller,
                        tmax = tmax_labeller)
  )
save_multiple(psd, "figures/two_loc_step_diff_sd_by_tmax")


## classification figure

two_loc_step_diff_classified <- readRDS("two_loc_step_diff_classified.rds")


loc1_change <- paste(two_loc_step_diff_classified$rt_ref_l1,
                     two_loc_step_diff_classified$rt_post_step_l1,
                     sep = " -> ")
loc2_change <- paste(two_loc_step_diff_classified$rt_ref_l2,
                     two_loc_step_diff_classified$rt_post_step_l2,
                     sep = " -> ")


two_loc_step_diff_classified$rt_change <- paste("Location 1: ", loc1_change,
                                         ", Location 2: ", loc2_change)
two_loc_step_diff_classified$true_eps <- factor(
  two_loc_step_diff_classified$true_eps, levels = eps_vals, ordered = TRUE
)
idx1 <- which(two_loc_step_diff_classified$true_label == "No transmission advantage" &
                two_loc_step_diff_classified$est_class == "Unclear")
idx2 <- which(two_loc_step_diff_classified$true_label == two_loc_step_diff_classified$est_class)
x <- two_loc_step_diff_classified[c(idx1, idx2), ]

p <- classification_fig(x)
save_multiple(
  p, "figures/two_loc_step_diff_classification"
)

if (! is.null(dev.list())) dev.off()
