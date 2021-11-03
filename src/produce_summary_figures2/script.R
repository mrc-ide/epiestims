## orderly::orderly_develop_start()
source("R/fig_utils.R")
dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
si_mu_ref <- 5.4
si_std_ref <- 1.5
##cols <- brewer.pal(8, "Set1")
values <- c(
  "Moderate" = "#ffa500",
  "Low" = "#005900",
  "High" = "#b20000",
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

error_summary <- map(infiles, readRDS)
## Appends multiplier label
error_summary <- affix_label(error_summary)
error_summary[["same_si"]] <- error_summary[["vary_si"]]
error_summary <- map2(
  error_summary, ms_vars, function(x, vars) {
    x[x$label %in% vars, ]
  }
)

## Label scenarios as low, central, high
error_summary <- map(error_summary, scenario_type_labeller)
error_summary <- map(
  error_summary, function(x) {
    x$metric <- "Bias"
    x
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
sd_summary[["same_si"]] <- sd_summary[["vary_si"]]
sd_summary <- map2(
  sd_summary, ms_vars, function(x, vars) {
    x[x$label %in% vars, ]
  }
)

sd_summary <- map(sd_summary, scenario_type_labeller)
sd_summary <- map(
  sd_summary, function(x) {
    x$metric <- "Uncertainty"
    x
})




classified <- readRDS("classification_by_scenario.rds")
classified <- affix_label(classified)
classified[["same_si"]] <- classified[["vary_si"]]
classified <- map2(
  classified, ms_vars, function(x, vars) {
    x[x$label %in% vars, ]
  }
)

classified <- map(classified, scenario_type_labeller)
classified <- map(
  classified, function(x) {
    x$metric <- "Classification"
    x
  }
)


## Same figures for SD
infiles <- list(
  vary_si = "vary_si_eps_summary_by_all_vars.rds",
  wrong_si = "wrong_si_eps_summary_by_all_vars.rds",
  vary_cv = "vary_cv_eps_summary_by_all_vars.rds",
  wrong_cv = "wrong_cv_eps_summary_by_all_vars.rds",
  vary_offs = "vary_offs_eps_summary_by_all_vars.rds",
  underrep = "underrep_eps_summary_by_all_vars.rds"
)

eps_summary <- map(infiles, readRDS)
eps_summary <- affix_label(eps_summary)
eps_summary[["same_si"]] <- eps_summary[["vary_si"]]
eps_summary <- map2(
  eps_summary, ms_vars, function(x, vars) {
    x[x$label %in% vars, ]
  }
)

eps_summary <- map(eps_summary, scenario_type_labeller)
eps_summary <- map(
  eps_summary, function(x) {
    x$metric <- "Coverage probability"
    x
})

scenarios <- names(eps_summary)
names(scenarios) <- scenarios
together <- map(
  scenarios, function(x) {
    x1 <- error_summary[[x]]
    x2 <- sd_summary[[x]]
    x3 <- eps_summary[[x]]
    x3 <- rename(
      x3, low = lower, med = pt_est, high = upper
    )
    x3 <- x3[, colnames(x2)]
    x4 <- classified[[x]]
    idx1 <- which(x4$true_label == "No transmission advantage" &
                  x4$est_class == "Unclear")
    idx2 <- which(x4$true_label == x4$est_class)
    x4 <- x4[c(idx1, idx2), ]
    x4 <- rename(
      x4, med = PointEst, low = Lower, high = Upper
    )
    x4 <- x4[, colnames(x2)]
    rbind(x1, x2, x3, x4)
  }
)


iwalk(
  together, function(x, scenario) {
    x <- split(x, list(x$tmax, x$rt_ref))
    iwalk(x, function(y, index) {
      y$scenario_type <- factor(
        y$scenario_type,
        levels = c("Baseline", "Low", "Moderate", "High"),
        ordered = TRUE
      )
      y$true_eps <- factor(
        y$true_eps, levels = unique(y$true_eps)
      )
      y$metric <- factor(
        y$metric, levels = c("Bias", "Uncertainty",
                             "Coverage probability",
                             "Classification"),
        ordered = TRUE
      )
      ## Coverage probability to go from 0 to 1
      dummy <- data.frame(
        metric = "Coverage probability",
        true_eps = levels(y$true_eps),
        low = 0, high = 1
      )
      dummy2 <- data.frame(metric = "Bias", y = 0)
      dummy$metric <- factor(
        dummy$metric, levels = levels(y$metric)
      )
      dummy2$metric <- factor(
        dummy2$metric, levels = levels(y$metric)
      )
      p <- ggplot(y) +
        geom_point(
          aes(true_eps, med, col = scenario_type),
          position = position_dodge(width = dodge_width),
          size = 1.4
        ) +
        geom_linerange(
          aes(true_eps, ymin = low, ymax = high,
              col = scenario_type),
          position = position_dodge(width = dodge_width),
          size = 1.1
        ) +
        geom_blank(
          data = dummy, aes(y = low)
        ) +
        geom_blank(
          data = dummy, aes(y = high)
        ) +
        geom_hline(
          data = dummy2, aes(yintercept = y),
          linetype = "dashed"
        ) +
        facet_wrap(~metric, scales = "free_y", ncol = 2) +
        theme_manuscript() +
        scale_color_manual(
          values = values,
          breaks = c("Low", "Moderate", "High")
        ) + labs(color = "Scenario Type") +
        xlab("Transmission Advantage") +
        ylab("")

      save_multiple(
        p, glue("figures/{scenario}_{index}")
      )
    }
    )
  }
)

