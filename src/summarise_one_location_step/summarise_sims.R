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

saveRDS(incid_summary, "incid_summary.rds")


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


eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), round, 3
)
eps_summary_df$true_eps <- round(eps_summary_df$true_eps, 3)

saveRDS(eps_summary_df, "eps_summary_df.rds")


## Summarise by parameters that vary
# ## 1. by tmax
# by_tmax <- split(eps_summary_df, eps_summary_df$tmax) %>%
#   map_dfr(
#     function(x) summarise_sims(na.omit(x)), .id = "tmax"
#   )
## 2. by epsilon
by_eps <- split(eps_summary_df, eps_summary_df$true_eps) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "true_eps"
  )

saveRDS(by_eps, "eps_summary_by_eps.rds")

## 6. by rt_ref and epsilon

## alternative summary function to generate 50% coverage
## df is grouped df, output of group_by
summarise_sims <- function(df) {
  summarise(
    df,
    n = sum(true_eps >= `2.5%` & true_eps <= `97.5%`),
    total = n(),
    pt_est = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 1],
    lower = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 2],
    upper = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 3],
    n50 = sum(true_eps >= `25%` & true_eps <= `75%`),
    ##total50 = n(),
    pt_est50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 1],
    lower50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 2],
    upper50 = Hmisc::binconf(x = n50, n = total, alpha = 0.05)[1, 3]
  )
}

by_eps_with_rt_change <- split(
  eps_summary_df, list(eps_summary_df$true_eps, eps_summary_df$rt_ref),
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

saveRDS(by_eps_with_rt_change, "eps_summary_by_eps_with_rt_change.rds")

## 7. by all variables (rt_ref, rt_post_step, tmax)
by_all_vars <-  split(
  eps_summary_df,
  list(eps_summary_df$rt_ref, eps_summary_df$rt_post_step, eps_summary_df$tmax,
       eps_summary_df$true_eps),
  sep = "_"
) %>%
  map_dfr(
    function(x) summarise_sims(na.omit(x)), .id = "var"
  )

by_all_vars <- tidyr::separate(by_all_vars, col = "var",
                               into = c("rt_ref", "rt_post_step", "tmax", "true_eps"),
                               sep = "_")

saveRDS(by_all_vars, "eps_summary_by_all_vars.rds")


# Error summary

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

saveRDS(eps_err_summary, "eps_err_summary.rds")


x <- as.list(sim_params)
x <- append(x, list(summary = eps_err_summary))


eps_err_summary_df <- pmap_dfr(
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
eps_err_summary_df <- na.omit(eps_err_summary_df)

saveRDS(eps_err_summary_df, "err_summary_df.rds")


# summarise the mean here (mu) but can also choose to summarise median (`50%`)
x <- group_by(eps_err_summary_df, rt_ref, rt_post_step, true_eps, tmax) %>%
    summarise(
      low = mean(mu) - sd(mu), med = mean(mu),
      high = mean(mu) + sd(mu)
    )

saveRDS(x, "err_summary_by_all_vars.rds")


x <- group_by(eps_err_summary_df, rt_ref, rt_post_step, true_eps, tmax) %>%
  summarise(
    low = mean(sd) - sd(sd), med = mean(sd),
    high = mean(sd) + sd(sd)
  )

saveRDS(x, "err_sd_summary_by_all_vars.rds")