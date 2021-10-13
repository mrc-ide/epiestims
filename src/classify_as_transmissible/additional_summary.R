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
  x <- split(x, x$confidence) %>%
    map_dfr(
      function(y) {
        tabyl(y, true_label, est_class, tmax) %>%
          adorn_percentages("row") %>%
          bind_rows(.id = "tmax")
      },.id = "confidence"
    )

  gather(
    x, classification, val, `Unclear`:`Variant more transmissible`
  )
}

## More coarse summary
## x is a subset of the classified data.frame
## holding a set of variables constant.
summary_other <- function(x) {
  tabyl(x, true_label, est_class) %>%
    adorn_percentages("row")
}

round_to <- 3
#################################################
#################################################
####### VARY SI #######
#################################################
#################################################
vary_si_eps_summary <- readRDS("vary_si_eps_summary_df.rds")
vary_si_eps_summary <- vary_si_eps_summary[vary_si_eps_summary$confidence != "Low", ]
vary_si_eps_summary <- mutate_at(
  vary_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_si_eps_summary$true_eps <- round(
  vary_si_eps_summary$true_eps, round_to
)
## What should be the correct classification,
## based on true advantage.
vary_si_eps_summary$true_label <- true_class(vary_si_eps_summary)
classified <- classify_epsilon(vary_si_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "vary_si_classified.rds")

## Summary across SI Mean
y <- split(classified, classified$si_mu_variant) %>%
  map_dfr(summary_other, .id = "Variant SI Mean")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_si_by_simu_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_si_by_tmax_classification.tex"
)


#################################################
#################################################
####### VARY SI #######
#################################################
#################################################
wrong_si_eps_summary <- readRDS("wrong_si_eps_summary_df.rds")
wrong_si_eps_summary <- wrong_si_eps_summary[wrong_si_eps_summary$confidence != "Low", ]
wrong_si_eps_summary <- mutate_at(
  wrong_si_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
wrong_si_eps_summary$true_eps <- round(
  wrong_si_eps_summary$true_eps, round_to
)
## What should be the correct classification,
## based on true advantage.
wrong_si_eps_summary$true_label <- true_class(wrong_si_eps_summary)
classified <- classify_epsilon(wrong_si_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "wrong_si_classified.rds")

## Summary across SI Mean
y <- split(classified, classified$si_mu_variant) %>%
  map_dfr(summary_other, .id = "Variant SI Mean")

cat(
  stargazer(y, summary = FALSE),
  file = "wrong_si_by_simu_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "wrong_si_by_tmax_classification.tex"
)

#################################################
#################################################
####### VARY CV #######
#################################################
#################################################
vary_cv_eps_summary <- readRDS("vary_cv_eps_summary_df.rds")
vary_cv_eps_summary <- vary_cv_eps_summary[vary_cv_eps_summary$confidence != "Low", ]
vary_cv_eps_summary <- mutate_at(
  vary_cv_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_cv_eps_summary$true_eps <- round(
  vary_cv_eps_summary$true_eps, round_to
)
vary_cv_eps_summary$true_label <- true_class(vary_cv_eps_summary)
classified <- classify_epsilon(vary_cv_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "vary_cv_classified.rds")

## Summary across SI CV
y <- split(classified, classified$si_cv_variant) %>%
  map_dfr(summary_other, .id = "Variant SI CV")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_cv_by_sicv_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_cv_by_tmax_classification.tex"
)



#################################################
#################################################
####### WRONG CV #######
#################################################
#################################################
wrong_cv_eps_summary <- readRDS("wrong_cv_eps_summary_df.rds")
wrong_cv_eps_summary <- wrong_cv_eps_summary[wrong_cv_eps_summary$confidence != "Low", ]
wrong_cv_eps_summary <- mutate_at(
  wrong_cv_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
wrong_cv_eps_summary$true_eps <- round(
  wrong_cv_eps_summary$true_eps, round_to
)
wrong_cv_eps_summary$true_label <- true_class(wrong_cv_eps_summary)
classified <- classify_epsilon(wrong_cv_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "wrong_cv_classified.rds")

## Summary across SI CV
y <- split(classified, classified$si_cv_variant) %>%
  map_dfr(summary_other, .id = "Variant SI CV")

cat(
  stargazer(y, summary = FALSE),
  file = "wrong_cv_by_sicv_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "wrong_cv_by_tmax_classification.tex"
)

#############################################
########### vary offspring
#############################################
vary_offs_eps_summary <- readRDS("vary_offs_eps_summary_df.rds")
vary_offs_eps_summary <- vary_offs_eps_summary[vary_offs_eps_summary$confidence != "Low", ]
vary_offs_eps_summary <- mutate_at(
  vary_offs_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
vary_offs_eps_summary$true_eps <- round(
  vary_offs_eps_summary$true_eps, round_to
)
vary_offs_eps_summary$true_label <- true_class(vary_offs_eps_summary)
classified <- classify_epsilon(vary_offs_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "vary_offs_classified.rds")

## Summary across SI CV
y <- split(classified, classified$kappa) %>%
  map_dfr(summary_other, .id = "Dispersion")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_offs_by_kappa_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "vary_offs_by_tmax_classification.tex"
)

#############################################
########### Under-reporting
#############################################
underrep_eps_summary <- readRDS("underrep_eps_summary_df.rds")
underrep_eps_summary <- underrep_eps_summary[underrep_eps_summary$confidence != "Low", ]
underrep_eps_summary <- mutate_at(
  underrep_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
underrep_eps_summary$true_eps <- round(
  underrep_eps_summary$true_eps, round_to
)
underrep_eps_summary$true_label <- true_class(underrep_eps_summary)
classified <- classify_epsilon(underrep_eps_summary)
tall <- summary_tmax_eps(classified)
saveRDS(tall, "underrep_classified.rds")

## Summary across SI CV
y <- split(classified, classified$p_report) %>%
  map_dfr(summary_other, .id = "Reporting probability")

cat(
  stargazer(y, summary = FALSE),
  file = "underrep_by_preport_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "underrep_by_tmax_classification.tex"
)
