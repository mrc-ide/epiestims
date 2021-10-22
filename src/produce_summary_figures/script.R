## orderly::orderly_develop_start()
## Aesthetics
source("R/fig_utils.R")
dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
cols <- brewer.pal(8, "Set1")
values <- c(
  "Central" = cols[5],
  "Low" = cols[3],
  "High" = cols[1],
  "Baseline" = "#222222"
)

## We have run more scenarios than we want to show in the
## manuscript. We have to filter out each separately
ms_vars <- list(
  vary_si = c("X 0.5", "X 1.5", "X 2"),
  wrong_si = c("X 0.5", "X 1.5", "X 2"),
  vary_cv = c("X 0.5", "X 1.5", "X 2"),
  wrong_cv = c("X 0.5", "X 1.5", "X 2"),
  vary_offs = c("0.1", "0.5", "1"),
  underrep = c("0.2", "0.5", "0.8"),
  same_si = "X 1"
)

scenario_names <- c(
  "same_si" = "(A) Baseline",
  "vary_offs" = "(B) With superspreading",
  "vary_si" = "(C) Different SI Mean",
  "wrong_si" = "(D) Misspecified SI Mean",
  "vary_cv" = "(E) Different SI CV",
  "wrong_cv" = "(F) Misspecified SI CV"
)

infiles <- list(
  vary_si = "vary_si_err_summary_by_all_vars.rds",
  wrong_si = "wrong_si_err_summary_by_all_vars.rds",
  vary_cv = "vary_cv_err_summary_by_all_vars.rds",
  wrong_cv = "wrong_cv_err_summary_by_all_vars.rds",
  vary_offs = "vary_offs_err_summary_by_all_vars.rds",
  underrep = "underrep_err_summary_by_all_vars.rds"
)

scenario_type_labeller <- function(x) {
  x$scenario_type <- case_when(
    x$label == "X 0.5" ~ "Low",
    x$label == "X 1.5" ~ "Central",
    x$label == "X 1" ~ "Baseline",
    x$label == "X 2" ~ "High",
    ## Vary offspring
    x$label == "0.1" ~ "High",
    x$label == "0.5" ~ "Central",
    x$label == "1" ~ "Low",
    ## Underreporting
    x$label == "0.2" ~ "High",
    x$label == "0.5" ~ "Central",
    x$label == "0.8" ~ "Low"
  )
  x
}
error_summary <- map(infiles, readRDS)
## Appends multiplier label
error_summary <- affix_label(error_summary)
## Label scenarios as low, central, high
error_summary <- map(error_summary, scenario_type_labeller)

split_df <- main_and_suppl(error_summary, ms_vars, ms_tmax)
## Split vary SI outputs into same SI and different SI
main_text_scenarios <- c("same_si", "vary_offs", "vary_si",
                         "wrong_si", "vary_cv", "wrong_cv")

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

error_summary$scenario_type <- factor(
  error_summary$scenario_type,
  levels = c("Baseline", "Low", "Central", "High"),
  ordered = TRUE
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
    aes(true_eps, med, col = scenario_type),
    position = position_dodge(width = dodge_width),
    size = 1.4
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high, col = scenario_type),
    position = position_dodge(width = dodge_width),
      size = 1.1
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
      breaks = c("Low", "Central", "High")
  ) + labs(color = "Scenario Type") +
  xlab("True Transmssion Advantage") +
  ylab("Bias") +
  theme_manuscript() +
  theme(legend.position = "bottom")

  save_multiple(p, glue("figures/main_text_fig_{index}"))
})

## Supplementary figures; error by tmax
iwalk(split_df, function(x, index) {
  p <- suppl_figure(x$suppl, index) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Error in estimating transmssion advantage")
  save_multiple(p, glue("figures/suppl_fig_{index}"))
})
## Same figures for SD
infiles <- list(
  vary_si = "vary_si_err_sd_summary_by_all_vars.rds",
  wrong_si = "wrong_si_err_sd_summary_by_all_vars.rds",
  vary_cv = "vary_cv_err_sd_summary_by_all_vars.rds",
  wrong_cv = "wrong_cv_err_sd_summary_by_all_vars.rds",
  vary_offs = "vary_offs_err_sd_summary_by_all_vars.rds",
  underrep = "underrep_err_sd_summary_by_all_vars.rds"
)

sd_summary <- map(infiles, readRDS)
sd_summary <- affix_label(sd_summary)
sd_summary <- map(sd_summary, scenario_type_labeller)
## give a fake ms_tmax so that everything goes to suppl
split_df <- main_and_suppl(sd_summary, ms_vars, ms_tmax = "60")
iwalk(split_df, function(x, index) {yes
  p <- suppl_figure(x$suppl, index) +
    ylab("Uncertainty")
  save_multiple(p, glue("figures/suppl_sd_fig_{index}"))
})
## Classification
classified <- readRDS("classification_by_scenario.rds")
classified <- split(classified, classified$scenario)
classified <- affix_label(classified)
plots <- iwalk(
  classified, function(x, scenario) {
    ggplot(x) +
      geom_point(aes(true_eps, n))
  }
)


if (! is.null(dev.list())) dev.off()
