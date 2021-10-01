## orderly::orderly_develop_start()
## Aesthetics
## df is a grouped dataframe with column med which is the
## median error
source("R/fig_utils.R")
summarise_median_err <- function(df, round_to = 3) {
  x <- summarise(
    df, median_low = quantile(med, 0.025),
    median_med = quantile(med, 0.5),
    median_high = quantile(med, 0.975)
  ) %>% ungroup()
  x <- mutate_if(x, is.numeric, round, round_to)
  x
}
## df is the output of summarise_median_err
format_median_err <- function(df) {
  df$formatted <-
  glue("{df$median_med}",
       " ({df$median_low}, {df$median_high})")
  df <- select(df, true_eps, tmax, formatted) %>%
    spread(key = tmax, value = formatted)
  df
}
## ms_vars is the set of variables we want to
## present in the main text. ms_tmax is the tmax
## value to use for results in main text
split_err_df <- function(x, ms_vars, ms_tmax = "50") {
  main <- x[x$label %in% ms_vars, ]
  main <- main[main$tmax == ms_tmax, ]

  suppl <- x[! x$label %in% ms_vars, ]
  suppl <- suppl[suppl$tmax != ms_tmax, ]
  list(main = main, suppl = suppl)
}

dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
round_to <- 3 ## Number of digits to round to

infiles <- list(
  vary_si = "vary_si_err_summary_by_all_vars.rds",
  wrong_si = "wrong_si_err_summary_by_all_vars.rds",
  vary_cv = "vary_cv_err_summary_by_all_vars.rds",
  wrong_cv = "wrong_cv_err_summary_by_all_vars.rds",
  vary_offs = "vary_offs_err_summary_by_all_vars.rds",
  underrep = "underrep_err_summary_by_all_vars.rds"
)

ms_vars <- list(
  vary_si = c("X 0.5", "X 1.5", "X 2"),
  wrong_si = c("X 0.5", "X 1.5", "X 2"),
  vary_cv = c("X 0.5", "X 1.5", "X 2"),
  wrong_cv = c("X 0.5", "X 1.5", "X 2"),
  vary_offs = c("0.1"),
  underrep = c("0.2", "0.5", "0.8"),
  same_si = "X 1"
)

error_summary <- map(infiles, readRDS)

error_summary[[1]]$label <- multiplier_label(
  error_summary[[1]]$si_mu_variant, si_mu_ref
)
error_summary[[2]]$label <- multiplier_label(
  error_summary[[2]]$si_mu_variant, si_mu_ref
)
error_summary[[3]]$label <- multiplier_label(
  error_summary[[3]]$si_cv_variant, si_std_ref / si_mu_ref
)
error_summary[[4]]$label <- multiplier_label(
  error_summary[[4]]$si_cv_variant, si_std_ref / si_mu_ref
)
error_summary[[5]]$label <- factor(error_summary[[5]]$kappa)
error_summary[[6]]$label <- factor(error_summary[[6]]$p_report)

## Split vary SI outputs into same SI and different SI
x <- error_summary[["vary_si"]]
same_si <- x[x$label == "X 1", ]
diff_si <- x[x$label != "X 1", ]
## Replace variable SI outputs in the list above
error_summary[["vary_si"]] <- diff_si
error_summary[["same_si"]] <- same_si

scenarios <- names(ms_vars)
names(scenarios) <- scenarios
split_df <- map(scenarios, function(sc) {
  split_err_df(error_summary[[sc]], ms_vars[[sc]])
})
main_text_scenarios <- c("same_si", "vary_offs", "vary_si",
                         "wrong_si", "vary_cv", "wrong_cv")
values <- c(
  "X 1" = "#222222", ## Same SI
  "0.1" = "#222222", ## Offspring
  ## Variable SI and CV.
  "X 0.5" = "#E69F00",
  "X 1.5" = "#56B4E9",
  "X 2" = "#009E73"
)
scenario_names <- c(
  "same_si" = "(A) Baseline",
  "vary_offs" = "(B) With superspreading",
  "vary_si" = "(C) Different SI Mean",
  "wrong_si" = "(D) Misspecified SI Mean",
  "vary_cv" = "(E) Different SI CV",
  "wrong_cv" = "(F) Misspecified SI CV"
)

main_text_df <- map_dfr(
  split_df[main_text_scenarios],
  ~ .[["main"]], .id = "scenario"
)

main_text_df$scenario <- factor(
  main_text_df$scenario,
  levels = main_text_scenarios,
  ordered = TRUE
)
main_text_df$true_eps <- factor(
  main_text_df$true_eps,
  levels = unique(main_text_df$true_eps)
)


main_text_df <- split(main_text_df, main_text_df$rt_ref)

iwalk(main_text_df, function(x, index) {
  ## Find the range for each scenario
  err_range <- group_by(x, scenario) %>%
    summarise(low = min(low), high = max(high)) %>%
    ungroup()
  ## Set the range to be same for all except
  ## wrong SI
  lowest <- min(err_range$low[err_range$scenario != "wrong_si"])
  highest <- max(err_range$high[err_range$scenario != "wrong_si"])
  err_range$low[err_range$scenario != "wrong_si"] <- floor(lowest)
  err_range$high[err_range$scenario != "wrong_si"] <- ceiling(highest)

  p <- ggplot(x) +
  geom_point(
    aes(true_eps, med, col = label),
      position = position_dodge(width = dodge_width),
      size = 1.4
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high, col = label),
      position = position_dodge(width = dodge_width),
      size = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(
    data = err_range, aes(y = low)
  ) +
  geom_blank(
    data = err_range, aes(y = high)
  ) +
  facet_wrap(
    ~scenario, scales = "free_y", ncol = 2,
    labeller = labeller(scenario = scenario_names)
  ) +
  scale_color_manual(
    values = values,
    breaks = c("X 0.5", "X 1.5", "X 2")
  ) + labs(color = "SI Mean or CV Multiplier") +
  xlab("True Transmssion Advantage") +
  ylab("Error in estimating transmssion advantage") +
  theme_manuscript() +
  theme(legend.position = "bottom")

  save_multiple(p, glue("figures/main_text_fig_{index}"))
})



## Classification
vary_si_classified <- readRDS("vary_si_classified.rds")
##vary_si_classified <- vary_si_classified[vary_si_classified$confidence != "Low", ]
p <- classification_fig(vary_si_classified)
save_multiple(p, "figures/vary_si_classification")

wrong_si_classified <- readRDS("wrong_si_classified.rds")
p <- classification_fig(wrong_si_classified)
save_multiple(p, "figures/wrong_si_classification")

vary_offs_classified <- readRDS("vary_offs_classified.rds")
p <- classification_fig(vary_offs_classified)
save_multiple(p, "figures/vary_offs_classification")

vary_cv_classified <- readRDS("vary_cv_classified.rds")
vary_cv_classified$true_eps <- as.numeric(vary_cv_classified$true_eps)
p <- classification_fig(vary_cv_classified)
save_multiple(p, "figures/vary_cv_classification")


wrong_cv_classified <- readRDS("wrong_cv_classified.rds")
wrong_cv_classified$true_eps <- as.numeric(wrong_cv_classified$true_eps)
p <- classification_fig(wrong_cv_classified)
save_multiple(p, "figures/wrong_cv_classification")

underrep_classified <- readRDS("underrep_classified.rds")
p <- classification_fig(underrep_classified)
save_multiple(p, "figures/underrep_classification")



if (! is.null(dev.list())) dev.off()

