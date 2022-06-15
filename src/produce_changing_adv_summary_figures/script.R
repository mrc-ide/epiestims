## orderly::orderly_develop_start(use_draft = "newer")
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
#ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
round_to <- 3 ## Number of digits to round to


infiles <- list(
  "10" = "err_summary_10.rds",
  "7" = "err_summary_7.rds",
  "standard" = "err_summary_standard.rds"
)

error_summary <- map(infiles, readRDS)

error_summary <- readRDS("one_loc_changing_adv_err_summary_df.rds")
error_summary$metric <- "Bias"

sd_summary <- readRDS("one_loc_changing_adv_err_sd_summary_df.rds")
sd_summary$metric <- "Uncertainty"

classified <- readRDS("classified.rds")
classified$metric <- "Classification"

eps_summary <- readRDS("one_loc_changing_adv_eps_summary_by_all_vars.rds")
eps_summary$metric <- "Coverage probability"


x1 <- error_summary
x1$qntl <- 'fake'
x2 <- sd_summary
x2$qntl <- 'fake'
x3 <- eps_summary
## This now has 50% coverage probability as well
x31 <- select(x3, tmax:upper, label:metric)
x32 <- select(
  x3, tmax, n = n50, pt_est = pt_est50,
  lower = lower50, upper = upper50, label:metric
)
x31$qntl <- '95%'
x32$qntl <- '50%'
x3 <- rbind(x31, x32)
x3 <- rename(
  x3, low = lower, med = pt_est, high = upper
)
x3 <- x3[, colnames(x2)]
x4 <- classified[[x]]
x4$qntl <- 1
idx1 <- which(x4$true_label == "No transmission advantage" &
                x4$est_class == "Unclear")
idx2 <- which(x4$true_label == x4$est_class)
x4 <- x4[c(idx1, idx2), ]
x4 <- rename(
  x4, med = PointEst, low = Lower, high = Upper
)
x4 <- x4[, colnames(x2)]
rbind(x1, x2, x3, x4)



one_loc_step_err <- readRDS("one_loc_step_err_summary_by_all_vars.rds")
eps_vals <- unique(one_loc_step_err$true_eps)

one_loc_step_err$true_eps <- factor(
  one_loc_step_err$true_eps, levels = eps_vals, ordered = TRUE
)
one_loc_step_err$rt_change <- factor(paste(one_loc_step_err$rt_ref,
                                           one_loc_step_err$rt_post_step,
                                           sep = "_to_")
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
                                          sep = "_to_")
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
                                                  sep = "_to_")
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
                                                sep = "_to_")
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

panel_fig(together, "one_loc_step")

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
                                           sep = "_to_")
)

two_loc_step_err$label <- multiplier_label(
  two_loc_step_err$si_mu_variant, si_mu_ref
)
two_loc_step_err$metric <- "Bias" # variable for faceting plots

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
                                          sep = "_to_")
)

two_loc_step_sd$label <- multiplier_label(
  two_loc_step_sd$si_mu_variant, si_mu_ref
)
two_loc_step_sd$metric <- "Uncertainty"

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
                                                  sep = "_to_")
)
two_loc_step_classified$true_eps <- factor(
  two_loc_step_classified$true_eps, levels = eps_vals, ordered = TRUE
)
idx1 <- which(two_loc_step_classified$true_label == "No transmission advantage" &
                two_loc_step_classified$est_class == "Unclear")
idx2 <- which(two_loc_step_classified$true_label == two_loc_step_classified$est_class)
x <- two_loc_step_classified[c(idx1, idx2), ]
x$metric <- "Classification"
y <- split(x, x$rt_change)

iwalk(y, function(change, index) {
  
  p <- classification_fig(change)
  save_multiple(
    p, glue("figures/two_loc_step_classification_{index}")
  )
  
})

two_loc_step_classified <- rename(x,
                                  med = PointEst, low = Lower, high = Upper)


# Coverage probability
two_loc_step_coverage <- readRDS("two_loc_step_eps_summary_by_all_vars.rds")
two_loc_step_coverage <- filter(two_loc_step_coverage, pt_est != "NaN")
two_loc_step_coverage$rt_change <- factor(paste(two_loc_step_coverage$rt_ref_l1,
                                                two_loc_step_coverage$rt_post_step_l1,
                                                sep = "_to_")
)
two_loc_step_coverage$true_eps <- factor(
  two_loc_step_coverage$true_eps, levels = eps_vals, ordered = TRUE
)
two_loc_step_coverage$metric <- "Coverage probability"
two_loc_step_coverage <- rename(two_loc_step_coverage,
                                med = pt_est, low = lower, high = upper)

two_loc_step_coverage$rt_ref_l1 <- as.numeric(two_loc_step_coverage$rt_ref_l1)
two_loc_step_coverage$rt_post_step_l1 <- as.numeric(two_loc_step_coverage$rt_post_step_l1)
two_loc_step_coverage$rt_ref_l2 <- as.numeric(two_loc_step_coverage$rt_ref_l2)
two_loc_step_coverage$rt_post_step_l2 <- as.numeric(two_loc_step_coverage$rt_post_step_l2)
two_loc_step_coverage$step_time_l1 <- as.numeric(two_loc_step_coverage$step_time_l1)
two_loc_step_coverage$step_time_l2 <- as.numeric(two_loc_step_coverage$step_time_l2)

together <- bind_rows(two_loc_step_err, two_loc_step_sd,
                      two_loc_step_coverage, two_loc_step_classified)

together$metric <- factor(
  together$metric, levels = c("Bias", "Uncertainty",
                              "Coverage probability",
                              "Classification"),
  ordered = TRUE
)

panel_fig(together, "two_loc_step")

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
two_loc_step_diff_err$metric <- "Bias" # variable for faceting plots

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
two_loc_step_diff_sd$metric <- "Uncertainty"


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
x$metric <- "Classification"

p <- classification_fig(x)
save_multiple(
  p, "figures/two_loc_step_diff_classification"
)

two_loc_step_diff_classified <- rename(x,
                                       med = PointEst, low = Lower, high = Upper)

# Coverage probability
two_loc_step_diff_coverage <- readRDS("two_loc_step_diff_eps_summary_by_all_vars.rds")
two_loc_step_diff_coverage <- filter(two_loc_step_diff_coverage, pt_est != "NaN")
two_loc_step_diff_coverage$rt_change <- factor(paste(two_loc_step_diff_coverage$rt_ref_l1,
                                                     two_loc_step_diff_coverage$rt_post_step_l1,
                                                     sep = " to ")
)
two_loc_step_diff_coverage$true_eps <- factor(
  two_loc_step_diff_coverage$true_eps, levels = eps_vals, ordered = TRUE
)
two_loc_step_diff_coverage$metric <- "Coverage probability"
two_loc_step_diff_coverage <- rename(two_loc_step_diff_coverage,
                                     med = pt_est, low = lower, high = upper)

two_loc_step_diff_coverage$rt_ref_l1 <- as.numeric(two_loc_step_diff_coverage$rt_ref_l1)
two_loc_step_diff_coverage$rt_post_step_l1 <- as.numeric(two_loc_step_diff_coverage$rt_post_step_l1)
two_loc_step_diff_coverage$rt_ref_l2 <- as.numeric(two_loc_step_diff_coverage$rt_ref_l2)
two_loc_step_diff_coverage$rt_post_step_l2 <- as.numeric(two_loc_step_diff_coverage$rt_post_step_l2)
two_loc_step_diff_coverage$step_time_l1 <- as.numeric(two_loc_step_diff_coverage$step_time_l1)
two_loc_step_diff_coverage$step_time_l2 <- as.numeric(two_loc_step_diff_coverage$step_time_l2)

together <- bind_rows(two_loc_step_diff_err, two_loc_step_diff_sd,
                      two_loc_step_diff_coverage, two_loc_step_diff_classified)

together$metric <- factor(
  together$metric, levels = c("Bias", "Uncertainty",
                              "Coverage probability",
                              "Classification"),
  ordered = TRUE
)

together$rt_change <- "1.4to1.1,1.6to1.2"

panel_fig(together, "two_loc_step_diff")

if (! is.null(dev.list())) dev.off()
