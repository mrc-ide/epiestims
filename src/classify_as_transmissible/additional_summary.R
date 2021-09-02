## orderly::orderly_develop_start()
## df is a summary data.frame with columns
## 2.5% and 97.5%
classify_epsilon <- function(df) {
  df$bias <- case_when(
    df$`2.5%` > 1 ~ "low_greater_than_1",
    df$`97.5%` < 1 ~ "high_less_than_1",
    TRUE ~ "CrI_includes_1"
  )
  df
}

vary_si_eps_summary <- readRDS("vary_si_eps_summary_df.rds")
vary_si_eps_summary <- mutate_at(
  vary_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, 3
)
vary_si_eps_summary$true_eps <- round(
  vary_si_eps_summary$true_eps, 3
)
## Perhaps this is only important when 95% CrI
## does not contain the true value
missed_true <- which(
  vary_si_eps_summary$true_eps >= vary_si_eps_summary$`2.5%` &
  vary_si_eps_summary$true_eps <= vary_si_eps_summary$`97.5%`
)
## True value not in 95% CrI in 64986 rows
## True value not in 95% CrI proportion 0.062
missed_true <- vary_si_eps_summary[-missed_true, ]
missed_true <- classify_epsilon(missed_true)

x <- tabyl(missed_true, true_eps, bias, tmax) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "tmax")

tall <- gather(x, label, val, -true_eps, -tmax)

saveRDS(tall, "vary_si_classified.rds")

## Summary across tmax
y <- tabyl(missed_true, tmax, bias) %>%
  adorn_percentages("row")
cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_si_by_tmax_classification.tex"
)

## Summary across true epsilon values
y <- tabyl(missed_true, true_eps, bias) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_si_by_eps_classification.tex"
)
