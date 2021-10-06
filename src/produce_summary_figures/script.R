## orderly::orderly_develop_start()
## Aesthetics
## df is a grouped dataframe with column med which is the
## median error
source("R/fig_utils.R")
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
## We have run more scenarios than we want to show in the
## manuscript. We have to filter out each separately
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

## Split into main and supplementart dfs
split_df <- list(
  same_si = list(
    main = filter(
      error_summary[["vary_si"]], label %in% ms_vars[["same_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["vary_si"]], label %in% ms_vars[["same_si"]],
      tmax != ms_tmax
    )
  ),
  vary_offs = list(
    main = filter(
      error_summary[["vary_offs"]], label %in% ms_vars[["vary_offs"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["vary_offs"]], label %in% ms_vars[["vary_offs"]],
      tmax != ms_tmax
    )
  ),
  vary_si = list(
    main = filter(
      error_summary[["vary_si"]], label %in% ms_vars[["vary_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["vary_si"]], label %in% ms_vars[["vary_si"]],
      tmax != ms_tmax
    )
  ),
  wrong_si = list(
    main = filter(
      error_summary[["wrong_si"]], label %in% ms_vars[["wrong_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["wrong_si"]], label %in% ms_vars[["wrong_si"]],
      tmax != ms_tmax
    )
  ),
  vary_cv = list(
    main = filter(
      error_summary[["vary_cv"]], label %in% ms_vars[["vary_cv"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["vary_cv"]], label %in% ms_vars[["vary_cv"]],
      tmax != ms_tmax
    )
  ),
  wrong_cv = list(
    main = filter(
      error_summary[["wrong_cv"]], label %in% ms_vars[["wrong_cv"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      error_summary[["wrong_cv"]], label %in% ms_vars[["wrong_cv"]],
      tmax != ms_tmax
    )
  ),
  underrep = list(
    main = NA, ## All of this goes in the supplementary
    suppl = filter(
      error_summary[["underrep"]], label %in% ms_vars[["underrep"]]
    )
  )
)
## Split vary SI outputs into same SI and different SI
main_text_scenarios <- c("same_si", "vary_offs", "vary_si",
                         "wrong_si", "vary_cv", "wrong_cv")
values <- c(
  "X 1" = "#222222", ## Same SI
  "0.1" = "#222222", ## Offspring
  "0.5" = "#D55E00",
  "1.0" = "#CC79A7",
  ## Variable SI and CV.
  "X 0.5" = "#E69F00",
  "X 1.5" = "#56B4E9",
  "X 2" = "#009E73",
  ## underreporting
  "0.2" = "#0072B2",
  "0.5" = "#D55E00",
  "0.8" = "#CC79A7"
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

## Supplementary figures; error by tmax
iwalk(split_df, function(x, index) {
  y <- x$suppl
  limits <- intersect(y$label, names(values))
  y$true_eps <- factor(
    y$true_eps,
    levels = unique(y$true_eps)
  )
  p <- ggplot(y) +
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
  facet_grid(
    tmax~rt_ref,
    labeller = labeller(tmax = tmax_labeller,
                        rt_ref = rt_labeller)
  ) +
  scale_color_manual(values = values, breaks = limits) +
  xlab("True Transmssion Advantage") +
  ylab("Error in estimating transmssion advantage") +
  theme_manuscript() +
    theme(legend.position = "bottom")
  if (index == "same_si") {
    p <- p + theme(legend.position = "none")
  } else if (index %in% c("vary_si", "vary_cv",
                          "wrong_si", "wrong_cv")) {
    p <- p + labs(color = "SI Mean or CV Multiplier")
  } else if (index == "vary_offs") {
    p <- p + labs(color = "Overdispersion")
  } else {
    ## That leaves under-reporting
    p <- p + labs(color = "Reporting probability")
  }
  save_multiple(p, glue("figures/suppl_fig_{index}"))
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

