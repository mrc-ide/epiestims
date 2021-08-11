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

## One of the problematic ones
##bad <- output[[105]][[5]][[98]]
##bad <- output[[1]][[5]][[68]]
x <- as.list(sim_params)
x <- append(x, list(summary = eps_summary))

eps_summary_df <- pmap_dfr(
  x, function(rt_ref, rt_post_step, step_time, epsilon,
              si_mu_variant, si_std_variant, summary) {
    summary$rt_ref <- rt_ref
    summary$rt_post_step <- rt_post_step
    summary$step_time <- step_time
    summary$true_eps <- epsilon
    summary$si_mu_variant <- si_mu_variant
    summary
  }
)

## idx <- which(
##   sim_params$rt_ref == sim_params$rt_ref[1] &
##   sim_params$si_mu_variant == sim_params$si_mu_variant[130]
## )
## x <- as.list(sim_params[idx, ])
## x <- append(x, list(summary = eps_summary[idx]))

## bad_summary_df <- pmap_dfr(
##   x, function(rt_ref, epsilon, si_mu_variant, si_std_variant, summary) {
##     summary$rt_ref <- rt_ref
##     summary$true_eps <- epsilon
##     summary$si_mu_variant <- si_mu_variant
##     summary
##   }
## )



eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), round, 3
)
eps_summary_df$true_eps <- round(eps_summary_df$true_eps, 3)

## bad_summary_df <- mutate_at(
##   bad_summary_df, vars(`2.5%`:`sd`), round, 3
## )
##bad_summary_df$true_eps <- round(bad_summary_df$true_eps, 3)

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
## 6. by rt_ref and rt_post_step
by_rt_change <- split(
  eps_summary_df, list(eps_summary_df$rt_ref, eps_summary_df$rt_post_step),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )
## 7. by all variables (rt_ref, rt_post_step, tmax)
by_all_vars <-  split(
  eps_summary_df,
  list(eps_summary_df$rt_ref, eps_summary_df$rt_post_step, eps_summary_df$tmax),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_all_vars <- tidyr::separate(by_all_vars, col = "var",
                               into = c("rt_ref", "rt_post_step", "tmax"), sep = "_")
# by_all_vars$si_mu_variant <- as.numeric(by_all_vars$si_mu_variant)
# by_all_vars$si_label <- glue(
#   "X {round(by_all_vars$si_mu_variant / si_mu_ref, 1)}"
# )
by_all_vars$rt_ref <- as.factor(by_all_vars$rt_ref)
by_all_vars$rt_post_step <- as.factor(by_all_vars$rt_post_step)

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

