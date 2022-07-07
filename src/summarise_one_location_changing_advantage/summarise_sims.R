## orderly::orderly_develop_start(parameters = list(estimation_window = 10))
## see comments in previous tasks
source("R/utils.R")
dir.create("figures")
unzip("output.zip")
sim_params <- readRDS("param_grid.rds")
sim_params$estimation_window <- estimation_window # add the window length for estimation epsilon
incid <- readRDS("incid.rds")


## Summarise simulated incidence, summarise
## epsilon and epsilon error.
tmax_all <- seq(10, 100, by = 10)
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
## for each row of params (nb here we only use one param 
## set so list has one element).
## Each element is a list of length 10, corresponding
## to the 10 tmax values we use (tmin is set to 19).
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
  x, function(rt_ref, epsilon_init, epsilon_final, epsilon_change,
              si_mu_variant, si_std_variant, estimation_window, summary) {
    summary$rt_ref <- rt_ref
    summary$epsilon_init <- epsilon_init
    summary$epsilon_final <- epsilon_final
    summary$epsilon_change <- epsilon_change
    summary$si_mu_variant <- si_mu_variant
    summary$estimation_window <- estimation_window
    summary
  }
)

eps_summary_df <- mutate_at(
  eps_summary_df, vars(`2.5%`:`sd`), round, 3
)

# add true_eps values
ndays <- 100
epsilon <- c(rep(sim_params$epsilon_init, sim_params$epsilon_change),
             seq(sim_params$epsilon_init, sim_params$epsilon_final, length.out = 30),
             rep(sim_params$epsilon_final, ndays - 30 - sim_params$epsilon_change))

epsilon <- data.frame(time = 1:100, true_eps = epsilon) %>% 
  filter(time %in% seq(19, 89, by = 10))

epsilon$time <- as.character(epsilon$time - 9)

eps_summary_df <- left_join(eps_summary_df, epsilon, by = c("tmax" = "time"))


eps_summary_df$true_eps <- round(eps_summary_df$true_eps, 3)

saveRDS(eps_summary_df, "eps_summary_df.rds")

# have the summary dataframe now! continue from here


# Error summary

eps_err_summary <- map(
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
    map2_dfr(
      res, names(res), function(res_tmax, t) {
        true_eps <- epsilon$true_eps[epsilon$time == t]
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon_error(res_sim[[1]], true_eps)
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
  x, function(rt_ref, epsilon_init, epsilon_final, epsilon_change,
              si_mu_variant, si_std_variant, estimation_window, summary) {
    summary$rt_ref <- rt_ref
    summary$epsilon_init <- epsilon_init
    summary$epsilon_final <- epsilon_final
    summary$epsilon_change <- epsilon_change
    summary$si_mu_variant <- si_mu_variant
    summary$estimation_window <- estimation_window
    summary
  }
)


eps_err_summary_df <- na.omit(eps_err_summary_df)

saveRDS(eps_err_summary_df, "eps_err_summary_df.rds")


# summarise the mean here (mu) but can also choose to summarise median (`50%`)
x <- group_by(eps_err_summary_df, tmax) %>%
  summarise(
    low = mean(mu) - sd(mu), med = mean(mu),
    high = mean(mu) + sd(mu)
  )

saveRDS(x, "err_summary_by_all_vars.rds")


x <- group_by(eps_err_summary_df, tmax) %>%
  summarise(
    low = mean(sd) - sd(sd), med = mean(sd),
    high = mean(sd) + sd(sd)
  )

saveRDS(x, "err_sd_summary_by_all_vars.rds")

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

by_all_vars <-  eps_summary_df %>%
  mutate(true_eps = replace_na(true_eps, 1.5)) %>% 
  group_by(tmax) %>%
  summarise_sims

saveRDS(by_all_vars, "eps_summary_by_all_vars.rds")