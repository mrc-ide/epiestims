## This is common to all summary tasks
dir.create("figures")
zip::unzip("output.zip")
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
    out <- map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon(res_sim[[1]])
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
    out <- estimate_uncertain(out)
    out
  }
)

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


for (i in 1:nrow(sim_params)) {
  eps_summary[[i]] <- cbind(
    eps_summary[[i]],
    sim_params[i, ]
  )
}

eps_summary_df <- do.call(
  what = 'rbind', args = eps_summary
)
eps_summary_df <- rename(eps_summary_df, true_eps = "epsilon")
message("eps_summary_df OK")
saveRDS(eps_summary_df, "eps_summary_df.rds")


for (i in 1:nrow(sim_params)) {
  eps_err_summary[[i]] <- cbind(
    eps_err_summary[[i]],
    sim_params[i, ]
  )
}

eps_err_summary_df <- do.call(what = 'rbind', args = eps_err_summary)
eps_err_summary_df <- rename(eps_err_summary_df, true_eps = "epsilon")
message("eps_err_summary_df OK")
saveRDS(eps_err_summary_df, "eps_err_summary_df.rds")

message("var = ", var)
x <- group_by(eps_err_summary_df, rt_ref, .data[[var]], true_eps, tmax) %>%
  summarise(
    low = mean(mu) - sd(mu), med = mean(mu),
    high = mean(mu) + sd(mu)
  ) %>% ungroup()

message("Summarising eps_err_summary_df by", var, "OK")
saveRDS(x, "err_summary_by_all_vars.rds")


x <- group_by(eps_err_summary_df, rt_ref, .data[[var]], true_eps, tmax) %>%
  summarise(
    low = mean(sd) - sd(sd), med = mean(sd),
    high = mean(sd) + sd(sd)
  ) %>% ungroup()

message("Summarising eps_err_summary_df SD by", var, "OK")
saveRDS(x, "err_sd_summary_by_all_vars.rds")


by_all_vars <-  group_by(
  eps_summary_df, rt_ref, .data[[var]], true_eps, tmax
) %>% summarise_sims

message("Summarising eps_summary_df by", var, "OK")
saveRDS(by_all_vars, "eps_summary_by_all_vars.rds")
