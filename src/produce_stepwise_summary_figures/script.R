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
# one_loc_step_coverage <- rename(one_loc_step_coverage,
                                # med = pt_est, low = lower, high = upper)

one_loc_step_coverage$rt_ref <- as.numeric(one_loc_step_coverage$rt_ref)
one_loc_step_coverage$rt_post_step <- as.numeric(one_loc_step_coverage$rt_post_step)

together <- list(error_summary = one_loc_step_err, sd_summary = one_loc_step_sd,
                 eps_summary = one_loc_step_coverage, classified = one_loc_step_classified)

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

together <- list(error_summary = two_loc_step_err, sd_summary = two_loc_step_sd,
                 eps_summary = two_loc_step_coverage, classified = two_loc_step_classified)


x1 <- together[[1]]
x1$qntl <- 'fake'
x2 <- together[[2]]
x2$qntl <- 'fake'
x3 <- together[[3]]
## This now has 50% coverage probability as well
x31 <- select(x3, rt_ref_l1:high, rt_change:metric)
x32 <- select(
  x3, rt_ref_l1:tmax, n = n50, total = total, med = pt_est50,
  low = lower50, high = upper50, rt_change:metric
)
x31$qntl <- '95%'
x32$qntl <- '50%'
x3 <- rbind(x31, x32)
x3$si_mu_variant <- 5.4
x3$label <- "X 1"
x3 <- x3[, colnames(x2)]
x4 <- together[[4]]
x4$qntl <- '1'
idx1 <- which(x4$true_label == "No transmission advantage" &
                x4$est_class == "Unclear")
idx2 <- which(x4$true_label == x4$est_class)
x4 <- x4[c(idx1, idx2), ]
x4$label <- "X 1"
x4 <- x4[, colnames(x2)]
out <- rbind(x1, x2, x3, x4)

out$metric <- factor(
  out$metric, levels = c("Bias", "Uncertainty",
                         "Coverage probability",
                         "Classification"),
  ordered = TRUE
)

## Construct dummy data.frame to control facet scales
## Coverage probability to go from 0 to 1
joined_data <- out
min_bias <- -1.5
max_bias <- 1.5
min_sd <- -0.1
max_sd <- 0.75
dummy <- data.frame(
  metric = c("Bias", "Coverage probability", "Uncertainty", "Classification"),
  ##true_eps = levels(y$true_eps),
  low = c(min_bias, 0, min_sd, 0),
  high = c(max_bias, 1, max_sd, 1)
)
dummy2 <- data.frame(
  metric = c("Bias", "Coverage probability",
             "Coverage probability"),
  y = c(0, 0.95, 0.5),
  qntl = c('fake', '95%', '50%')
)

dummy$metric <- factor(
  dummy$metric, levels = levels(joined_data$metric)
)
dummy2$metric <- factor(
  dummy2$metric, levels = levels(joined_data$metric)
)

## Different shapes for 95% and 50% coverage probability
joined_data$shape <- 19 ## everything is a circle
joined_data$shape[joined_data$qntl == "50%"] <- 18 ## Except 50% coverage probability


y <- split(joined_data, list(joined_data$tmax, joined_data$rt_change))
iwalk(y, function(z, index) {
  
  p <- ggplot(z) +
    geom_point(
      aes(true_eps, med, group = qntl, shape = shape), col = "black",
      position = position_dodge2(width = dodge_width),
      size = 1.2
    ) +
    geom_linerange(
      aes(true_eps, ymin = low, ymax = high, group = qntl),
      col = "black",
      position = position_dodge2(width = dodge_width)
    ) +
    geom_blank(
      data = dummy, aes(y = low)
    ) +
    geom_blank(
      data = dummy, aes(y = high)
    ) +
    geom_hline(
      data = dummy2, aes(yintercept = y, group = qntl),
      linetype = "dashed"
    ) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_shape_identity() +
    theme_manuscript() +
    xlab("Transmission Advantage") +
    ylab("")
  
  save_multiple(
    p, glue("figures/two_loc_step_panel_{index}")
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

together <- list(error_summary = two_loc_step_diff_err, sd_summary = two_loc_step_diff_sd,
                 eps_summary = two_loc_step_diff_coverage, classified = two_loc_step_diff_classified)

x1 <- together[[1]]
x1$qntl <- 'fake'
x2 <- together[[2]]
x2$qntl <- 'fake'
x3 <- together[[3]]
## This now has 50% coverage probability as well
x31 <- select(x3, rt_ref_l1:high, rt_change:metric)
x32 <- select(
  x3, rt_ref_l1:tmax, n = n50, total = total, med = pt_est50,
  low = lower50, high = upper50, rt_change:metric
)
x31$qntl <- '95%'
x32$qntl <- '50%'
x3 <- rbind(x31, x32)
x3$si_mu_variant <- 5.4
x3$label <- "X 1"
x3 <- x3[, colnames(x2)]
x4 <- together[[4]]
x4$qntl <- '1'
idx1 <- which(x4$true_label == "No transmission advantage" &
                x4$est_class == "Unclear")
idx2 <- which(x4$true_label == x4$est_class)
x4 <- x4[c(idx1, idx2), ]
x4$label <- "X 1"
x4 <- x4[, colnames(x2)]
out <- rbind(x1, x2, x3, x4)
out$rt_change <- "1.4to1.1,1.6to1.2"
out$metric <- factor(
  out$metric, levels = c("Bias", "Uncertainty",
                         "Coverage probability",
                         "Classification"),
  ordered = TRUE
)

## Construct dummy data.frame to control facet scales
## Coverage probability to go from 0 to 1
joined_data <- out
min_bias <- -1.5
max_bias <- 1.5
min_sd <- -0.1
max_sd <- 0.75
dummy <- data.frame(
  metric = c("Bias", "Coverage probability", "Uncertainty", "Classification"),
  ##true_eps = levels(y$true_eps),
  low = c(min_bias, 0, min_sd, 0),
  high = c(max_bias, 1, max_sd, 1)
)
dummy2 <- data.frame(
  metric = c("Bias", "Coverage probability",
             "Coverage probability"),
  y = c(0, 0.95, 0.5),
  qntl = c('fake', '95%', '50%')
)

dummy$metric <- factor(
  dummy$metric, levels = levels(joined_data$metric)
)
dummy2$metric <- factor(
  dummy2$metric, levels = levels(joined_data$metric)
)

## Different shapes for 95% and 50% coverage probability
joined_data$shape <- 19 ## everything is a circle
joined_data$shape[joined_data$qntl == "50%"] <- 18 ## Except 50% coverage probability


y <- split(joined_data, list(joined_data$tmax, joined_data$rt_change))
iwalk(y, function(z, index) {
  
  p <- ggplot(z) +
    geom_point(
      aes(true_eps, med, group = qntl, shape = shape), col = "black",
      position = position_dodge2(width = dodge_width),
      size = 1.2
    ) +
    geom_linerange(
      aes(true_eps, ymin = low, ymax = high, group = qntl),
      col = "black",
      position = position_dodge2(width = dodge_width)
    ) +
    geom_blank(
      data = dummy, aes(y = low)
    ) +
    geom_blank(
      data = dummy, aes(y = high)
    ) +
    geom_hline(
      data = dummy2, aes(yintercept = y, group = qntl),
      linetype = "dashed"
    ) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_shape_identity() +
    theme_manuscript() +
    xlab("Transmission Advantage") +
    ylab("")
  
  save_multiple(
    p, glue("figures/two_loc_step_diff_panel_{index}")
  )
})


if (! is.null(dev.list())) dev.off()
