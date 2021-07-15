## orderly::orderly_develop_start(use_draft = "newer")
## see comments in previous tasks
source("R/utils.R")
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

## Get 1 summary row for each row in the parameter grid
eps_summary_df <- pmap_dfr(
  list(summary = eps_summary, true_epsilon = params$epsilon),
  function(summary, true_epsilon) {
    summary <- mutate_at(
      summary, vars(`2.5%`:`sd`), round, digits = 2
    )
    true_epsilon <- round(true_epsilon, 3)
    out <- group_by(summary, tmax) %>%
      summarise(
        n = sum(true_epsilon >=`2.5%` & true_epsilon <= `97.5%`),
        total = n()
      ) %>% ungroup()
    pivot_wider(
      out, names_from = "tmax", values_from = c("n", "total")
    )
  }
)

eps_summary_df <- cbind(params, eps_summary_df)
## Now tall for plotting
eps_summary_tall <- pivot_longer(
  eps_summary_df, cols = n_20:n_60,
  names_to = c("param", "tmax"),
  names_sep = "_"
)
## total_20, total_30 etc should all be the same
## and equal to nsims.
ci <- pmap(
  eps_summary_tall[, c("value", "total_20")],
  function(value, total_20) {
    Hmisc::binconf(x = value, n = total_20, alpha = 0.05)
  }
)
eps_summary_tall$pt_est <- map_dbl(ci, ~ .[[1]])
eps_summary_tall$lower <- map_dbl(ci, ~ .[[2]])
eps_summary_tall$upper <- map_dbl(ci, ~ .[[3]])

## columns of interest
eps_summary_tall <- select(
  eps_summary_tall, rt_ref:si_std_variant, tmax:upper
)

eps_summary_tall$rt_ref <- factor(eps_summary_tall$rt_ref)
eps_summary_tall$si_mu_variant <- factor(eps_summary_tall$si_mu_variant)

ggplot(eps_summary_tall) +
  geom_point(aes(tmax, pt_est)) +
  geom_linerange(aes(x = tmax, ymin = lower, ymax = upper)) +
  facet_grid(si_mu_variant ~ epsilon)

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



saveRDS(eps_err_summary_df, "eps_err_summary.rds")
saveRDS(eps_summary_df, "eps_summary_df.rds")
saveRDS(incid_summary, "incid_summary.rds")

