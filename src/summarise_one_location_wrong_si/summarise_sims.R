## orderly::orderly_develop_start()
## see comments in previous tasks
prefix <- "wrong_si"
source("R/utils.R")
source("R/summarise_sims.R")
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
## To include in manuscript
by_all_vars$rt_ref <- factor(by_all_vars$rt_ref)
by_all_vars$true_eps <- factor(by_all_vars$true_eps, levels = eps_vals, ordered = TRUE)

ms_tmax <- "50"
ms_df <- by_all_vars[by_all_vars$tmax == ms_tmax, ]
si_df <- by_all_vars[by_all_vars$tmax != ms_tmax, ]
p1 <- true_epsilon_vs_95CrI(ms_df)
ggsave(glue("figures/{prefix}_by_true_eps_tmax_{ms_tmax}.png"), p1)

p2 <- true_epsilon_vs_95CrI(si_df) +
  facet_wrap(~tmax)

ggsave(glue("figures/{prefix}_by_true_eps_tmax_si.png"), p2)



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

x <- group_by(eps_err_summary_df, rt_ref, true_eps, tmax) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  )
x <- ungroup(x)
x$rt_ref <- as.factor(x$rt_ref)
x$true_eps <- factor(x$true_eps, levels = eps_vals, ordered = TRUE)



p1 <- true_epsilon_vs_error(x) +
  facet_wrap(~tmax)

ggsave(glue("figures/{prefix}_by_true_eps_error.png"), p1)


saveRDS(eps_err_summary, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")
