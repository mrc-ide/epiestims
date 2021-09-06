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
round_to <- 3
vary_si_eps_summary <- readRDS("vary_si_eps_summary_df.rds")
vary_si_eps_summary <- mutate_at(
  vary_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_si_eps_summary$true_eps <- round(
  vary_si_eps_summary$true_eps, round_to
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


vary_si_eps_summary <- readRDS("vary_si_eps_summary_df.rds")
vary_si_eps_summary <- mutate_at(
  vary_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_si_eps_summary$true_eps <- round(
  vary_si_eps_summary$true_eps, round_to
)
## ## Perhaps this is only important when 95% CrI
## ## does not contain the true value
## missed_true <- which(
##   vary_si_eps_summary$true_eps >= vary_si_eps_summary$`2.5%` &
##   vary_si_eps_summary$true_eps <= vary_si_eps_summary$`97.5%`
## )
## ## True value not in 95% CrI in 64986 rows
## ## True value not in 95% CrI proportion 0.062
## missed_true <- vary_si_eps_summary[-missed_true, ]
missed_true <- classify_epsilon(vary_si_eps_summary)

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



vary_si_eps_summary <- readRDS("vary_si_eps_summary_df.rds")
vary_si_eps_summary <- mutate_at(
  vary_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_si_eps_summary$true_eps <- round(
  vary_si_eps_summary$true_eps, round_to
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

#################################
#################################
####### VARY CV #######
#################################
#################################
vary_cv_eps_summary <- readRDS("vary_cv_eps_summary_df.rds")
vary_cv_eps_summary <- mutate_at(
  vary_cv_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_cv_eps_summary$true_eps <- round(
  vary_cv_eps_summary$true_eps, round_to
)
missed_true <- which(
  vary_cv_eps_summary$true_eps >= vary_cv_eps_summary$`2.5%` &
  vary_cv_eps_summary$true_eps <= vary_cv_eps_summary$`97.5%`
)
## True value not in 95% CrI in 3751 rows
## True value not in 95% CrI proportion 0.057
missed_true <- vary_cv_eps_summary[-missed_true, ]
missed_true <- classify_epsilon(missed_true)

x <- tabyl(missed_true, true_eps, bias, tmax) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "tmax")

tall <- gather(x, label, val, -true_eps, -tmax)

saveRDS(tall, "vary_cv_classified.rds")

## Summary across tmax
y <- tabyl(missed_true, tmax, bias) %>%
  adorn_percentages("row")
cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_cv_by_tmax_classification.tex"
)

## Summary across true epsilon values
y <- tabyl(missed_true, true_eps, bias) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_cv_by_eps_classification.tex"
)
#############################################
########### vary offspring
vary_offs_eps_summary <- readRDS("vary_offs_eps_summary_df.rds")
vary_offs_eps_summary <- mutate_at(
  vary_offs_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_offs_eps_summary$true_eps <- round(
  vary_offs_eps_summary$true_eps, round_to
)
missed_true <- which(
  vary_offs_eps_summary$true_eps >= vary_offs_eps_summary$`2.5%` &
  vary_offs_eps_summary$true_eps <= vary_offs_eps_summary$`97.5%`
)
## True value not in 95% CrI in 9792 rows
## True value not in 95% CrI proportion 0.251
missed_true <- vary_offs_eps_summary[-missed_true, ]
missed_true <- classify_epsilon(missed_true)

x <- tabyl(missed_true, true_eps, bias, tmax) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "tmax")

tall <- gather(x, label, val, -true_eps, -tmax)

saveRDS(tall, "vary_offs_classified.rds")

## Summary across tmax
y <- tabyl(missed_true, tmax, bias) %>%
  adorn_percentages("row")
cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_offs_by_tmax_classification.tex"
)

## Summary across true epsilon values
y <- tabyl(missed_true, true_eps, bias) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_offs_by_eps_classification.tex"
)
