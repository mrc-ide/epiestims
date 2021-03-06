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
#ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
round_to <- 3 ## Number of digits to round to

values <- c(
  "10 days" = "#ffa500",
  "All data" = "#b20000",
  "7 days" = "#005900")


infiles <- list(
  "10" = "err_summary_10.rds",
  "7" = "err_summary_7.rds",
  "standard" = "err_summary_standard.rds"
)

error_summary <- map(infiles, readRDS)
error_summary <- map2(
  error_summary, names(error_summary), function(x, y) {
    x$window <- y
    x$metric <- "Bias"
    x
  })


## Same figures for SD
infiles <- list(
  "10" = "err_sd_summary_10.rds",
  "7" = "err_sd_summary_7.rds",
  "standard" = "err_sd_summary_standard.rds"
)

sd_summary <- map(infiles, readRDS)
sd_summary <- map2(
  sd_summary, names(sd_summary), function(x, y) {
    x$window <- y
    x$metric <- "Uncertainty"
    x
  })

##Same for classification
infiles <- list(
  "10" = "classified_10.rds",
  "7" = "classified_7.rds",
  "standard" = "classified_standard.rds"
)

classified <- map(infiles, readRDS)
classified <- flatten(classified)
names(classified) <- names(infiles)
classified <- map2(
  classified, names(classified), function(x, y) {
    x$window <- y
    x$metric <- "Classification"
    x
  })

## And finally for coverage probability
infiles <- list(
  "10" = "eps_summary_by_all_vars_10.rds",
  "7" = "eps_summary_by_all_vars_7.rds",
  "standard" = "eps_summary_by_all_vars_standard.rds"
)

eps_summary <- map(infiles, readRDS)
eps_summary <- map2(
  eps_summary, names(eps_summary), function(x, y) {
    x$window <- y
    x$metric <- "Coverage probability"
    x
  })


scenarios <- names(eps_summary)
names(scenarios) <- scenarios

## Append a fake column "qntl" to all metrics to
## make plotting easier.
together <- map(
  scenarios, function(x) {
    x1 <- error_summary[[x]]
    x1$qntl <- 'fake'
    x2 <- sd_summary[[x]]
    x2$qntl <- 'fake'
    x3 <- eps_summary[[x]]
    ## This now has 50% coverage probability as well
    x31 <- select(x3, tmax:upper, window:metric)
    x32 <- select(
      x3, tmax, n = n50, total = total, pt_est = pt_est50,
      lower = lower50, upper = upper50, window:metric
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
    # browser()
    x4 <- x4[c(idx1, idx2), ]
    x4 <- rename(
      x4, med = PointEst, low = Lower, high = Upper
    )
    x4 <- x4[, colnames(x2)]
    rbind(x1, x2, x3, x4)
  }
)

## Produce panel figure

baseline <- bind_rows(together)

baseline$metric <- factor(
  baseline$metric, levels = c("Bias", "Uncertainty",
                              "Coverage probability",
                              "Classification"),
  ordered = TRUE
)

baseline$window <- factor(
  baseline$window,
  levels = c("standard", "10", "7"),
  labels = c("All data", "10 days", "7 days"),
  ordered = TRUE
)
breaks <- c("All data", "10 days", "7 days")
baseline <- arrange(baseline, window)

baseline$tmax <- as.numeric(baseline$tmax)
baseline <- filter(baseline, tmax < 90)
baseline$tmax <- factor(
  baseline$tmax, levels = unique(baseline$tmax)
)


## Construct dummy data.frame to control facet scales
## Coverage probability to go from 0 to 1
min_bias <- -1.5
max_bias <- 1.5
min_sd <- -0.25
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
  dummy$metric, levels = levels(baseline$metric)
)
dummy2$metric <- factor(
  dummy2$metric, levels = levels(baseline$metric)
)
## Different shapes for 95% and 50% coverage probability
baseline$shape <- 19 ## everything is a circle
baseline$shape[baseline$qntl == "50%"] <- 18 ## Except 50% coverage probability

p <- ggplot(baseline) +
  geom_point(
    aes(tmax, med, col = window, group = qntl, shape = shape),
    size = 1.5, position = position_dodge2(width = dodge_width)
  ) +
  geom_linerange(
    aes(tmax, ymin = low, ymax = high,
        col = window, group = qntl), position = position_dodge2(width = dodge_width)
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
  scale_color_manual(
    values = values,
    breaks = breaks
  ) +
  xlab("Time (days)") +
  ylab("")
p
p <- p + labs(color = "Window length")
p

save_multiple(p, "figures/one_loc_changing_adv")


## Create figure showing how the estimated epsilon varies over the scenario

infiles <- list(
  "10" = "eps_summary_df_10.rds",
  "7" = "eps_summary_df_7.rds",
  "standard" = "eps_summary_df_standard.rds"
)

eps_summary_df <- map(infiles, readRDS)
eps_summary_df <- map_dfr(
  eps_summary_df, function(x) {
    x$estimation_window <- as.character(x$estimation_window)
    x
  })

eps_summary_df$tmax <- as.numeric(eps_summary_df$tmax)
eps_summary_df <- arrange(eps_summary_df, tmax) %>% 
  filter(tmax < 90)

levels <- c("standard", "10", "7")
labels <- c("All data", "10 days", "7 days")
eps_summary_df$estimation_window <- factor(eps_summary_df$estimation_window,
                                           levels = levels, labels = labels)

# summarise across all sims for each tmax and window value
eps_summary_df <- group_by(eps_summary_df, estimation_window, tmax, true_eps) %>%
  summarise(
    low = mean(mu) - sd(mu), med = mean(mu),
    high = mean(mu) + sd(mu)
  )

# recreate the epsilon values used in the simulations
sim_params <- expand.grid(
  rt_ref = 1.1,
  epsilon_init = 1.1,
  epsilon_final = 1.5,
  epsilon_change = 30,
  si_mu_variant = 1 * si_mu_ref,
  si_std_variant = si_std_ref
)
ndays <- 100
epsilon_all <- c(rep(sim_params$epsilon_init, sim_params$epsilon_change),
                 seq(sim_params$epsilon_init, sim_params$epsilon_final, length.out = 30),
                 rep(sim_params$epsilon_final, ndays - 30 - sim_params$epsilon_change))

epsilon_all <- data.frame(time = c(1:9, (19:100) - 9),
                          epsilon = c(rep(1.1, 9), epsilon_all[19:100]))

# create plot
p <- ggplot(eps_summary_df) +
  geom_line(data = epsilon_all,
            aes(x = time, y = epsilon), linetype = "dashed") +
  geom_point(aes(x = tmax, y = med, colour = estimation_window),
             size = 1.5,
             position = position_dodge(width = 3)) +
  geom_linerange(aes(x = tmax, ymax = high, ymin = low, colour = estimation_window),
                 size = 1,
                 position = position_dodge(width = 3)) +
  geom_point(aes(x = tmax, y = true_eps),
             size = 1.5, shape = 0) +
  scale_y_continuous(limits = c(0, 1.8)) +
  scale_color_manual(name = "Window length",
                     values = values,
                     breaks = c("All data", "10 days", "7 days")) +
  xlab("Time (days)") +
  ylab("Effective transmission advantage") +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

p

save_multiple(p, "figures/one_loc_changing_adv_epsilon_est")

########################
## Create figure showing how the estimated epsilon varies over the scenario
## Use the aggregated posterior distribution

infiles <- list(
  "10" = "posterior_eps_summary_10.rds",
  "7" = "posterior_eps_summary_7.rds",
  "standard" = "posterior_eps_summary_standard.rds"
)

post_summary_df <- map(infiles, readRDS)
post_summary_df <- bind_rows(post_summary_df, .id = "estimation_window")

post_summary_df$tmax <- as.numeric(post_summary_df$tmax)
post_summary_df <- arrange(post_summary_df, tmax) %>% 
  filter(tmax < 90)

levels <- c("standard", "10", "7")
labels <- c("All data", "10 days", "7 days")
post_summary_df$estimation_window <- factor(post_summary_df$estimation_window,
                                           levels = levels, labels = labels)

# summarise across all sims for each tmax and window value
# # post_summary_df <- group_by(post_summary_df, estimation_window, tmax) %>%
#   summarise(
#     low = mu - sd, med = mu,
#     high = mu + sd
#   )

# recreate the epsilon values used in the simulations
sim_params <- expand.grid(
  rt_ref = 1.1,
  epsilon_init = 1.1,
  epsilon_final = 1.5,
  epsilon_change = 30,
  si_mu_variant = 1 * si_mu_ref,
  si_std_variant = si_std_ref
)
ndays <- 100
epsilon_all <- c(rep(sim_params$epsilon_init, sim_params$epsilon_change),
                 seq(sim_params$epsilon_init, sim_params$epsilon_final, length.out = 30),
                 rep(sim_params$epsilon_final, ndays - 30 - sim_params$epsilon_change))

epsilon_all <- data.frame(time = c(1:9, (19:100) - 9),
                          epsilon = c(rep(1.1, 9), epsilon_all[19:100]))

# create plot
p <- ggplot(post_summary_df) +
  geom_line(data = epsilon_all,
            aes(x = time, y = epsilon), linetype = "dashed") +
  geom_point(aes(x = tmax, y = median, colour = estimation_window),
             size = 1.5,
             position = position_dodge(width = 3)) +
  geom_linerange(aes(x = tmax, ymax = `97.5%`, ymin = `2.5%`, colour = estimation_window),
                 size = 1,
                 position = position_dodge(width = 3)) +
  geom_point(data = eps_summary_df, aes(x = tmax, y = true_eps),
             size = 1.5, shape = 0) +
  scale_y_continuous(limits = c(0, 2.4)) +
  scale_color_manual(name = "Window length",
                     values = values,
                     breaks = c("All data", "10 days", "7 days")) +
  xlab("Time (days)") +
  ylab("Effective transmission advantage") +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

p

p

save_multiple(p, "figures/one_loc_changing_adv_posterior_summary")

if (! is.null(dev.list())) dev.off()
