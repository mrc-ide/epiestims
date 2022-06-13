## orderly::orderly_develop_start()
## df is a summary data.frame with columns
## 2.5% and 97.5%
classify_epsilon <- function(df) {
  df$est_class <- case_when(
    df$`2.5%` > 1 ~ "Variant more transmissible",
    df$`97.5%` < 1 ~ "Variant less transmissible",
    TRUE ~ "Unclear"
  )
  df
}

true_class <- function(df) {
  case_when(
    ## In our case this is when eps is actually 1
    df$true_eps <= 1 ~ "No transmission advantage",
    ## All other scenarios are "More transmissible"
    ## so the CrI should comfortably exclude 1.
    TRUE ~ "Variant more transmissible"
  )
}

## Summarise classification performance as a
## function of tmax and true advantage.
summary_tmax_eps <- function(x) {
  ## This method has the advantage
  ## that numbers across the three
  ## estimated class will add to the
  ## total number of simulations
  out <- group_by(x, across(!sim & !est_class)) %>%
    summarise(
      Unclear = sum(est_class == 'Unclear'),
      `Variant more transmissible` = sum(est_class == "Variant more transmissible"),
      `Variant less transmissible` = sum(est_class == "Variant less transmissible")
    ) %>% ungroup()
  ## This will of course be the number of sims
  out$nsims <- out$Unclear +
    out$`Variant more transmissible` +
    out$`Variant less transmissible`

  out <- gather(
    out, est_class, n, Unclear:`Variant less transmissible`
  )

  ci_95 <- binconf(out$n, out$nsims) %>%
    tidy()

  cbind(out, ci_95[['x']])
}

## More coarse summary
## x is a subset of the classified data.frame
## holding a set of variables constant.
summary_other <- function(x) {
  tabyl(x, true_label, est_class) %>%
    adorn_percentages("row")
}


infiles <- list(
  vary_si = "vary_si_eps_summary_df.rds",
  wrong_si = "wrong_si_eps_summary_df.rds",
  vary_cv = "vary_cv_eps_summary_df.rds",
  wrong_cv = "wrong_cv_eps_summary_df.rds",
  vary_offs = "vary_offs_eps_summary_df.rds",
  underrep = "underrep_eps_summary_df.rds"
)

eps_summary <- map(infiles, readRDS)

classified <- map(
  eps_summary, function(x) {
    x$true_label <- true_class(x)
    classified <- classify_epsilon(x)
    x <- select(
      classified, tmax, sim, rt_ref:est_class
    )
    summary_tmax_eps(x)
  }
)

saveRDS(
  classified, "classification_by_scenario.rds"
)

classified <- bind_rows(classified, .id = "scenario")
write_csv(classified, "classification_by_scenario.csv")

eps1 <- classified[classified$true_label == "No transmission advantage", ]
## Distinguish baseline scenario
eps1$scenario[eps1$si_mu_variant == 5.4 & eps1$scenario == "vary_si"] <- "baseline"
x <- count(eps1, scenario, tmax, est_class, wt = n)
x <- spread(x, est_class, n)
write_csv(x, "classification_by_scenario_by_tamx_eps1.csv")

x <- count(eps1, scenario, est_class, wt = n)
write_csv(
  x, "classification_by_scenario_eps1.csv"
)
