## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
prefix <- "vary_si"
source("R/utils.R")
source("R/summarise_sims.R")
dir.create("figures")
unzip("output.zip")
x <- as.list(sim_params)
x <- append(x, list(summary = eps_summary))

eps_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, si_mu_variant, si_std_variant, summary) {
    summary$rt_ref <- rt_ref
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
## 1. by tmax
by_tmax <- split(eps_summary_df, eps_summary_df$tmax) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "tmax"
  )
## 2. by epsilon
by_eps <- split(eps_summary_df, eps_summary_df$true_eps) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "true_eps"
)
## 3. by si_mu
by_mu_var <- split(eps_summary_df, eps_summary_df$si_mu_variant) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "si_mu_variant"
  )
## 4. by rt_ref
by_rt_ref <- split(eps_summary_df, eps_summary_df$rt_ref) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "rt_ref"
  )
## 5. by rt_ref and si_mu
by_rt_and_mu <- split(
  eps_summary_df, list(eps_summary_df$rt_ref, eps_summary_df$si_mu_variant),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )
## 6. by all variables
by_all_vars <-  split(
  eps_summary_df,
  list(eps_summary_df$rt_ref, eps_summary_df$si_mu_variant, eps_summary_df$tmax),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_all_vars <- tidyr::separate(by_all_vars, col = "var", into = c("rt_ref", "si_mu_variant", "tmax"), sep = "_")
by_all_vars$si_mu_variant <- as.numeric(by_all_vars$si_mu_variant)
by_all_vars$si_label <- glue(
  "X {round(by_all_vars$si_mu_variant / si_mu_ref, 1)}"
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
  theme(
    legend.position = "top"
  )

ggsave(glue("figures/vary_si_prop_in_95.png"), p)

## Proportion in 95% CrI by tmax
by_all_vars <-  split(
  eps_summary_df,
  list(eps_summary_df$rt_ref, eps_summary_df$true_eps, eps_summary_df$tmax),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_all_vars <- tidyr::separate(by_all_vars, col = "var", into = c("rt_ref", "true_eps", "tmax"), sep = "_")
by_all_vars$rt_ref <- as.factor(by_all_vars$rt_ref)
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


## Error
x <- as.list(sim_params)
x <- append(x, list(summary = eps_err_summary))

eps_err_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, si_mu_variant, si_std_variant, summary) {
    summary$rt_ref <- rt_ref
    summary$true_eps <- epsilon
    summary$si_mu_variant <- si_mu_variant
    summary
  }
)
eps_err_summary_df <- na.omit(eps_err_summary_df)
x <- group_by(eps_err_summary_df, rt_ref, si_mu_variant, tmax) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  )

x$si_label <- glue(
  "X {round(x$si_mu_variant / si_mu_ref, 1)}"
)
x$rt_ref <- as.factor(x$rt_ref)

p <- ggplot(x) +
  geom_point(
    aes(tmax, med, col = factor(rt_ref)),
    position = position_dodge(width = 0.5)
  ) +
  geom_linerange(
    aes(tmax, ymin = low, ymax = high, col = factor(rt_ref)),
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~si_label) +
  theme_minimal() +
  labs(color = "Reference Rt") +
  theme(legend.position = "top")

ggsave(glue("figures/eps_error.png"), p)

saveRDS(eps_err_summary, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")

