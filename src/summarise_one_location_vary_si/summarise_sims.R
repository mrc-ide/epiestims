## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
source("R/utils.R")
dir.create("figures")
params <- readRDS("one_loc_vary_si_params.rds")
incid <- readRDS("one_loc_vary_si_incid.rds")
output <- readRDS("one_loc_vary_si.rds")

## Summarise simulated incidence, summarise
## epsilon and epsilon error.
tmax_all <- seq(20, 60, by = 10)
names(tmax_all) <- tmax_all

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
  output, function(res) {
    names(res) <- tmax_all
    map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon(res_sim)
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
  }
)

x <- as.list(params)
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
by_tmax <- summarise_sims(group_by(eps_summary_df, tmax))
## 2. by epsilon
by_eps <- summarise_sims(group_by(eps_summary_df, true_eps))
## 3. by si_mu
by_mu_var <- summarise_sims(group_by(eps_summary_df, si_mu_variant))
## 4. by rt_ref
by_rt_ref <- summarise_sims(group_by(eps_summary_df, rt_ref))
## 5. by rt_ref and si_mu
by_rt_and_mu <- summarise_sims(
  group_by(eps_summary_df, rt_ref, si_mu_variant)
)
## 6. by all variables
by_all_vars <- summarise_sims(
  group_by(eps_summary_df, rt_ref, si_mu_variant, tmax)
) %>% ungroup()

si_mu_ref <- 6.83
by_all_vars$si_label <- glue(
  "X {round(by_all_vars$si_mu_variant / si_mu_ref, 1)}"
)

p <- ggplot(by_all_vars) +
  geom_point(aes(tmax, pt_est)) +
  geom_linerange(aes(tmax, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Proportion in 95% CrI") +
  xlab("tmax") +
  ylim(0, 1)

##
p1 <- p +
  facet_grid_paginate(rt_ref ~ si_label, ncol = 2, nrow = 3, page = 1)
pages <- n_pages(p1)

walk(
  seq_len(pages), function(page) {
    p1 <- p +
      facet_grid_paginate(rt_ref ~ si_label, ncol = 2, nrow = 3, page = page)
    ggsave(glue("figures/vary_si_prop_in_95_{page}.png"), p1)
  }
)

eps_err_summary <- pmap(
  list(res = output, epsilon = params$epsilon),
  function(res, epsilon) {
    names(res) <- tmax_all
    map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon_error(res_sim, epsilon)
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
  }
)



saveRDS(eps_err_summary, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")

