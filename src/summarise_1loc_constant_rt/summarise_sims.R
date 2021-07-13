##orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
source("R/utils.R")
params <- readRDS("one_loc_constant_rt_params.rds")
incid <- readRDS("one_loc_constant_rt_incid.rds")
output <- readRDS("one_loc_constant_rt.rds")

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

eps_summary_df <- pmap_dfr(
  list(
    rt_ref = params$rt_ref, epsilon = params$epsilon,
    summary = eps_summary
  ),
  function(rt_ref, epsilon, summary) {
    summary$rt_ref <- rt_ref
    summary$true_epsilon <- epsilon
    summary
  }
)

eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), signif, digits = 3
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

eps_err_summary_df <- pmap_dfr(
  list(
    rt_ref = params$rt_ref, epsilon = params$epsilon,
    summary = eps_err_summary
  ),
  function(rt_ref, epsilon, summary) {
    summary$rt_ref <- rt_ref
    summary$true_epsilon <- epsilon
    summary
  }
)

eps_err_summary_df <- mutate_at(
  eps_err_summary_df, vars(`2.5%`:`sd`), signif, digits = 3
)

##
prop_in_95 <- group_by(eps_summary_df, tmax) %>%
  summarise(
    n = sum(true_epsilon >= `2.5%` & true_epsilon <= `97.5%`),
    total = n(),
    lower = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 2],
    upper = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 3]
  ) %>% ungroup()

saveRDS(eps_err_summary_df, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")
saveRDS(prop_in_95, "prop_in_95.rds")

