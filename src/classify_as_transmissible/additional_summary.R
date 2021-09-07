## orderly::orderly_develop_start()
## df is a summary data.frame with columns
## 2.5% and 97.5%
classify_epsilon <- function(df) {
  df$est_class <- case_when(
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
## What should be the correct classification,
## based on true advantage.
vary_si_eps_summary$true_class <- case_when(
  ## In our case this is when eps is actually 1
  vary_si_eps_summary$true_eps <= 1 ~ "CrI_includes_1",
  ## All other scenarios are "More transmissible"
  ## so the CrI should comfortably exclude 1.
  TRUE ~ "low_greater_than_1"
)

classified <- classify_epsilon(vary_si_eps_summary)
classified$classification <- case_when(
  classified$true_class == classified$est_class ~ "CORRECT",
  TRUE ~ "INCORRECT"
)

x <- tabyl(classified, tmax, true_eps, classification) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "classification")


tall <- gather(x, true_eps, val, -classification, -tmax)

saveRDS(tall, "vary_si_classified.rds")

## Summary across SI Mean
y <- tabyl(classified, si_mu_variant, classification) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_si_by_simu_classification.tex"
)

## Summary across tmax
y <- tabyl(classified, tmax, classification) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_si_by_tmax_classification.tex"
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

vary_cv_eps_summary$true_class <- case_when(
  ## In our case this is when eps is actually 1
  vary_cv_eps_summary$true_eps <= 1 ~ "CrI_includes_1",
  ## All other scenarios are "More transmissible"
  ## so the CrI should comfortably exclude 1.
  TRUE ~ "low_greater_than_1"
)

classified <- classify_epsilon(vary_cv_eps_summary)
classified$classification <- case_when(
  classified$true_class == classified$est_class ~ "CORRECT",
  TRUE ~ "INCORRECT"
)


x <- tabyl(classified, tmax, true_eps, classification) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "classification")


tall <- gather(x, true_eps, val, -classification, -tmax)
saveRDS(tall, "vary_cv_classified.rds")

## Summary across SI Mean
y <- tabyl(classified, si_cv_variant, classification) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_cv_by_sicv_classification.tex"
)

## Summary across tmax
y <- tabyl(classified, tmax, classification) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_cv_by_tmax_classification.tex"
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

vary_offs_eps_summary$true_class <- case_when(
  ## In our case this is when eps is actually 1
  vary_offs_eps_summary$true_eps <= 1 ~ "CrI_includes_1",
  ## All other scenarios are "More transmissible"
  ## so the CrI should comfortably exclude 1.
  TRUE ~ "low_greater_than_1"
)

classified <- classify_epsilon(vary_offs_eps_summary)
classified$classification <- case_when(
  classified$true_class == classified$est_class ~ "CORRECT",
  TRUE ~ "INCORRECT"
)


x <- tabyl(classified, tmax, true_eps, classification) %>%
  adorn_percentages("row") %>%
  bind_rows(.id = "classification")


tall <- gather(x, true_eps, val, -classification, -tmax)
saveRDS(tall, "vary_offs_classified.rds")

## Summary across SI Mean
y <- tabyl(classified, kappa, classification) %>%
  adorn_percentages("row")

cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_offs_by_kappa_classification.tex"
)

## Summary across tmax
y <- tabyl(classified, tmax, classification) %>%
  adorn_percentages("row")



cat(
  stargazer(y, summary = FALSE, row.names = FALSE),
  file = "vary_offs_by_tmax_classification.tex"
)

y <- tabyl(classified, tmax, kappa, classification) %>%
  adorn_percentages("row")

