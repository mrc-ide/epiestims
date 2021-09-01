## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
prefix <- "vary_offspring"
dir.create("figures")
source("R/utils.R")
source("R/summarise_sims.R")

x <- as.list(sim_params)
x <- append(x, list(summary = eps_summary))

eps_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, kappa, summary) {
    summary$rt_ref <- rt_ref
    summary$true_eps <- epsilon
    summary$kappa <- kappa
    summary
  }
)
saveRDS(eps_summary_df, "eps_summary_df.rds")

## 6. by all variables
by_all_vars <-  group_by(
  eps_summary_df,
  rt_ref, kappa, tmax, true_eps
) %>% summarise_sims

saveRDS(by_all_vars, "eps_summary_by_all_vars.rds")

## Error
x <- as.list(sim_params)
x <- append(x, list(summary = eps_err_summary))

eps_err_summary_df <- pmap_dfr(
  x, function(rt_ref, epsilon, kappa, summary) {
    summary$rt_ref <- rt_ref
    summary$true_eps <- epsilon
    summary$kappa <- kappa
    summary
  }
)

saveRDS(eps_err_summary_df, "eps_err_summary_df.rds")

x <- group_by(eps_err_summary_df, rt_ref, kappa, true_eps, tmax) %>%
  summarise(
    low = quantile(`50%`, 0.025), med = quantile(`50%`, 0.5),
    high = quantile(`50%`, 0.975)
  ) %>% ungroup()

saveRDS(x, "err_summary_by_all_vars.rds")

