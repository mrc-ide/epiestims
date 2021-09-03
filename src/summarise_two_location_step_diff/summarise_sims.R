## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
source("R/utils.R")
dir.create("figures")
unzip("output.zip")
sim_params <- readRDS("param_grid.rds")
incid <- readRDS("incid.rds")


## Summarise simulated incidence, summarise
## epsilon and epsilon error.
tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all
si_mu_ref <- 5.4
si_std_ref <- 1.5
si_distr_ref <- discr_si(0:30, mu = si_mu_ref, sigma = si_std_ref)
incid_summary <- map_depth(
  incid, 2, function(x) {
    map_dfr(
      tmax_all, function(tmax) {
        data.frame(
          ref_incid = sum(x[1:tmax, 1, 1]),
          var_incid = sum(x[1:tmax, 1, 2]),
          tmax = tmax
        )
      }
    )
  }
)

incid_summary <- map(
  incid_summary, ~ bind_rows(., .id = "sim")
)


## Structure of output: list with one element
## for each row of params.
## Each element is a list of length 5, corresponding
## to the 5 tmax values we use (tmin is set to 15).
## Each element at nested level 2 is a list of
## length 100, corresponding to the number of
## simulations that we run.
## Finally within that, we have a list of length
## 2, which is the output from estimate_joint.
eps_summary <- map(
  seq_len(nrow(sim_params)),
  function(index) {
    infile <- glue("outputs/estimate_joint_{index}.rds")
    if (!file.exists(infile)) {
      warning(infile, " not present")
      return(NULL)
    }
    res <- readRDS(infile)
    message("At ", infile)
    names(res) <- tmax_all
    map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon(res_sim[[1]])
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
  }
)

x <- as.list(sim_params)
x <- append(x, list(summary = eps_summary))

eps_summary_df <- pmap_dfr(
  x, function(rt_ref_l1, rt_post_step_l1, step_time_l1,
              rt_ref_l2, rt_post_step_l2, step_time_l2,
              epsilon, si_mu_variant, si_std_variant, summary) {
    summary$rt_ref_l1 <- rt_ref_l1
    summary$rt_post_step_l1 <- rt_post_step_l1
    summary$step_time_l1 <- step_time_l1
    summary$rt_ref_l2 <- rt_ref_l2
    summary$rt_post_step_l2 <- rt_post_step_l2
    summary$step_time_l2 <- step_time_l2
    summary$true_eps <- epsilon
    summary$si_mu_variant <- si_mu_variant
    summary
  }
)


eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), round, 3
)
eps_summary_df$true_eps <- round(eps_summary_df$true_eps, 3)


## Summarise by parameters that vary

## 1. by epsilon
by_eps <- split(eps_summary_df, eps_summary_df$true_eps) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "true_eps"
  )

## 2. by rt_ref and rt_post_step
by_eps_with_rt_change <- split(
  eps_summary_df, list(eps_summary_df$true_eps, eps_summary_df$rt_ref_l1),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_eps_with_rt_change <- tidyr::separate(by_eps_with_rt_change, col = "var",
                                         into = c("true_eps", "rt_ref"), sep = "_")

by_eps_with_rt_change$rt_ref <- as.factor(by_eps_with_rt_change$rt_ref)

eps_vals <- unique(by_eps_with_rt_change$true_eps)

by_eps_with_rt_change$true_eps <- factor(by_eps_with_rt_change$true_eps,
                                         levels = eps_vals, ordered = TRUE)

# Plot summarising over all tmax values

p <-
  ggplot(by_eps_with_rt_change) +
  geom_point(aes(true_eps, pt_est, colour = rt_ref),
             position = position_dodge(width = 0.5)) +
  geom_linerange(aes(true_eps, ymin = lower, ymax = upper, colour = rt_ref),
                 position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Proportion in 95% CrI") +
  xlab("True transmission advantage") +
  ylim(0, 1) +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(
    legend.position = "top"
  )

ggsave(glue("figures/two_loc_step_diff_prop_in_95.png"), p)


# Plot summarising only at tmax = 50

eps_summary_t50 <- eps_summary_df[eps_summary_df$tmax == "50", ]

by_eps_with_rt_change_t50 <- split(
  eps_summary_t50,
  list(eps_summary_t50$true_eps, eps_summary_t50$rt_ref_l1),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_eps_with_rt_change_t50 <- tidyr::separate(by_eps_with_rt_change_t50, col = "var",
                                         into = c("true_eps", "rt_ref"), sep = "_")

by_eps_with_rt_change_t50$rt_ref <- as.factor(by_eps_with_rt_change_t50$rt_ref)

eps_vals <- unique(by_eps_with_rt_change_t50$true_eps)

by_eps_with_rt_change_t50$true_eps <- factor(by_eps_with_rt_change_t50$true_eps,
                                         levels = eps_vals, ordered = TRUE)


p1 <-
  ggplot(by_eps_with_rt_change_t50) +
  geom_point(aes(true_eps, pt_est, colour = rt_ref),
             position = position_dodge(width = 0.5)) +
  geom_linerange(aes(true_eps, ymin = lower, ymax = upper, colour = rt_ref),
                 position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Proportion in 95% CrI") +
  xlab("True transmission advantage") +
  ylim(0, 1) +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(
    legend.position = "top"
  )

ggsave(glue("figures/two_loc_step_diff_prop_in_95_tmax50.png"), p1)




# Absolute errors

eps_err_summary <- map2(
  seq_len(nrow(sim_params)),
  sim_params$epsilon,
  function(index, epsilon) {
    infile <- glue("outputs/estimate_joint_{index}.rds")
    if (!file.exists(infile)) {
      warning(infile, " not present")
      return(NULL)
    }
    res <- readRDS(infile)
    message("At ", infile)
    names(res) <- tmax_all
    map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon_error(res_sim[[1]], epsilon)
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
  }
)

x <- as.list(sim_params)
x <- append(x, list(summary = eps_err_summary))

eps_err_summary_df <- pmap_dfr(
  x, function(rt_ref_l1, rt_post_step_l1, step_time_l1,
              rt_ref_l2, rt_post_step_l2, step_time_l2,
              epsilon, si_mu_variant, si_std_variant, summary) {
    summary$rt_ref_l1 <- rt_ref_l1
    summary$rt_post_step_l1 <- rt_post_step_l1
    summary$step_time_l1 <- step_time_l1
    summary$rt_ref_l2 <- rt_ref_l2
    summary$rt_post_step_l2 <- rt_post_step_l2
    summary$step_time_l2 <- step_time_l2
    summary$true_eps <- epsilon
    summary$si_mu_variant <- si_mu_variant
    summary
  }
)

# Errors across all tmax
eps_err_summary_df <- na.omit(eps_err_summary_df)
x <- group_by(eps_err_summary_df, rt_ref_l1, true_eps) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  )

x$rt_ref_l1 <- as.factor(x$rt_ref_l1)
eps_vals <- unique(x$true_eps)
x$true_eps <- factor(x$true_eps,
                     levels = eps_vals, ordered = TRUE)

p <-
  ggplot(x) +
  geom_point(
    aes(true_eps, med, colour = rt_ref_l1),
    position = position_dodge(width = 0.5)
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high, colour = rt_ref_l1),
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Estimated - True transmission advantage") +
  xlab("True transmission advantage") +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(legend.position = "top")

ggsave(glue("figures/eps_error.png"), p)


# Summarise absolute error at tmax = 50

# Errors across all tmax
eps_err_summary_df_t50 <- eps_err_summary_df[eps_summary_df$tmax == "50", ]

x <- group_by(eps_err_summary_df_t50, rt_ref_l1, true_eps) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  )

x$rt_ref_l1 <- as.factor(x$rt_ref_l1)
eps_vals <- unique(x$true_eps)
x$true_eps <- factor(x$true_eps,
                     levels = eps_vals, ordered = TRUE)

p2 <-
  ggplot(x) +
  geom_point(
    aes(true_eps, med, colour = rt_ref_l1),
    position = position_dodge(width = 0.5)
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high, colour = rt_ref_l1),
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Estimated - True transmission advantage") +
  xlab("True transmission advantage") +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(legend.position = "top")

ggsave(glue("figures/eps_error_tmax50.png"), p)




saveRDS(eps_err_summary, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")

