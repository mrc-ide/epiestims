## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
prefix <- "vary_cv"
source("R/utils.R")
source("R/summarise_sims.R")
x <- as.list(sim_params)
x <- append(x, list(summary = eps_summary))

eps_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, si_mu_variant, si_cv_variant, summary) {
    summary$rt_ref <- rt_ref
    summary$true_eps <- epsilon
    summary$si_cv_variant <- si_cv_variant
    summary
  }
)


eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), round, 3
)
eps_summary_df$true_eps <- round(eps_summary_df$true_eps, 3)


by_all_vars <-  group_by(
  eps_summary_df,
  rt_ref, si_cv_variant, tmax, true_eps
) %>% summarise_sims

by_all_vars$si_cv_variant <- as.numeric(by_all_vars$si_cv_variant)
si_cv_ref <- si_std_ref / si_mu_ref
by_all_vars$si_label <- glue(
  "X {round(by_all_vars$si_cv_variant / si_cv_ref, 1)}"
)
by_all_vars$rt_ref <- as.factor(by_all_vars$rt_ref)



p <- ggplot(by_all_vars) +
  geom_point(aes(tmax, pt_est, col = rt_ref)) +
  geom_linerange(aes(tmax, ymin = lower, ymax = upper, col = rt_ref)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Proportion in 95% CrI") +
  xlab("tmax") +
  ylim(0, 1) +
  facet_wrap(~ si_label, nrow = 3) +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(legend.position = "top")

ggsave(glue("figures/vary_cv_prop_in_95.png"), p)

## For MS. x-axis - true epsilon; y-axis: prop in 95% CrI
## tmax = fixed
by_all_vars <- group_by(
  eps_summary_df, rt_ref, tmax, true_eps
) %>% summarise_sims

by_all_vars <- ungroup(by_all_vars)
eps_vals <- unique(by_all_vars$true_eps)
split(by_all_vars, by_all_vars$tmax) %>%
  iwalk(
    function(x, tmax) {
      x$rt_ref <- factor(x$rt_ref)
      x$true_eps <- factor(x$true_eps, levels = eps_vals, ordered = TRUE)
      p <- ggplot(x) +
        geom_point(
          aes(true_eps, pt_est, col = rt_ref),
          position = position_dodge(width = 0.3)
        ) +
        geom_linerange(
          aes(true_eps, ymin = lower, ymax = upper, col = rt_ref),
          position = position_dodge(width = 0.3)
        ) +
        geom_hline(yintercept = 0.95, linetype = "dashed") +
        ylab("Proportion in 95% CrI") +
        xlab("True transmission advantage") +
        ylim(0, 1) +
        theme_minimal() +
        labs(color = "Reference Rt") +
        theme(
          legend.position = "top"
        )

      ggsave(glue("figures/{prefix}_by_true_eps_tmax_{tmax}.png"), p)
    }
  )


x <- as.list(sim_params)
x <- append(x, list(summary = eps_err_summary))

eps_err_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, si_mu_variant, si_cv_variant, summary) {
    summary$rt_ref <- rt_ref
    summary$true_eps <- epsilon
    summary$si_cv_variant <- si_cv_variant
    summary
  }
)
eps_err_summary_df <- na.omit(eps_err_summary_df)
x <- group_by(eps_err_summary_df, rt_ref, tmax, true_eps) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  )
x$rt_ref <- as.factor(x$rt_ref)

split(x, x$tmax) %>%
  iwalk(
    function(y, tmax) {
      y$rt_ref <- factor(y$rt_ref)
      y$true_eps <- factor(y$true_eps, levels = eps_vals, ordered = TRUE)
      p <- ggplot(y) +
        geom_point(
          aes(true_eps, med, col = rt_ref),
          position = position_dodge(width = 0.3)
        ) +
        geom_linerange(
          aes(true_eps, ymin = low, ymax = high, col = rt_ref),
          position = position_dodge(width = 0.3)
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Estimated - True transmission advantage ") +
        xlab("True transmission advantage") +
        theme_minimal() +
        labs(color = "Reference Rt") +
        theme(
          legend.position = "top"
        )

      ggsave(glue("figures/{prefix}_err_by_true_eps_tmax_{tmax}.png"), p)
    }
  )


saveRDS(eps_err_summary, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")

