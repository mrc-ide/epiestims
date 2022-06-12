## orderly::orderly_develop_start()
## x is a list of data.frames -
## either mean error or SD
affix_label <- function(x) {
 x[["vary_si"]]$label <- multiplier_label(
   x[["vary_si"]]$si_mu_variant, si_mu_ref
 )
 x[["wrong_si"]]$label <- multiplier_label(
   x[["wrong_si"]]$si_mu_variant, si_mu_ref
 )
 x[["vary_cv"]]$label <- multiplier_label(
   x[["vary_cv"]]$si_cv_variant, si_std_ref / si_mu_ref
 )
 x[["wrong_cv"]]$label <- multiplier_label(
   x[["wrong_cv"]]$si_cv_variant, si_std_ref / si_mu_ref
 )
 x[["vary_offs"]]$label <- factor( x[["vary_offs"]]$kappa)
 ## To distinguish it from vary offspring
 x[["underrep"]]$label <- paste0("p", x[["underrep"]]$p_report)
 x
}

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
  "Baseline" = "#222222",
  "0.5" = "#ffa500",
  "0.8" = "#005900",
  "0.2" = "#b20000"
)

## We have run more scenarios than we want to show in the
## manuscript. We have to filter out each separately
ms_vars <- list(
  vary_si = c("X 0.5", "X 1.5", "X 2"),
  wrong_si = c("X 0.5", "X 1.5", "X 2"),
  vary_cv = c("X 0.5", "X 1.5", "X 2"),
  wrong_cv = c("X 0.5", "X 1.5", "X 2"),
  vary_offs = c("0.1", "0.5", "1"),
  underrep = c("p0.2", "p0.5", "p0.8"),
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
## Append a fake column "qntl" to all metrics to
## make plotting easier.
together <- map(
  scenarios, function(x) {
    x1 <- error_summary[[x]]
    x1$qntl <- 'fake'
    x2 <- sd_summary[[x]]
    x2$qntl <- 'fake'
    x3 <- eps_summary[[x]]
    ## This now has 50% coverage probability as well
    x31 <- select(x3, rt_ref:upper, label:metric)
    x32 <- select(
      x3, rt_ref:tmax, n = n50, pt_est = pt_est50,
      lower = lower50, upper = upper50, label:metric
    )
    x31$qntl <- '95%'
    x32$qntl <- '50%'
    x3 <- rbind(x31, x32)
    x3 <- rename(
      x3, low = lower, med = pt_est, high = upper
    )
    x3 <- x3[, colnames(x2)]
    x4 <- classified[[x]]
    x4$qntl <- 1
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
      if (scenario == "underrep") {
        y$scenario_type <- factor(
          y$scenario_type,
          levels = c("0.2", "0.5", "0.8"),
          ordered = TRUE
        )
        breaks <- c("0.2", "0.5", "0.8")
      } else {
        y$scenario_type <- factor(
          y$scenario_type,
          levels = c("Baseline", "Low", "Moderate", "High"),
          ordered = TRUE
        )
        breaks <- c("Low", "Moderate", "High")
      }
      y$true_eps <- factor(
        y$true_eps, levels = unique(y$true_eps)
      )
      y$metric <- factor(
        y$metric, levels = c("Bias", "Uncertainty",
                             "Coverage probability",
                             "Classification"),
        ordered = TRUE
      )

      ## Construct dummy data.frame to control facet scales
      ## Coverage probability to go from 0 to 1
      min_bias <- ifelse(scenario == "wrong_si", -2.5, -1.5)
      max_bias <- ifelse(scenario == "wrong_si", 17.5, 1.5)
      min_sd <- -0.1
      max_sd <- ifelse(scenario == "wrong_si", 1, 0.75)
      dummy <- data.frame(
        metric = c("Bias", "Coverage probability", "Uncertainty"),
        ##true_eps = levels(y$true_eps),
        low = c(min_bias, 0, min_sd),
        high = c(max_bias, 1, max_sd)
      )
      dummy2 <- data.frame(
        metric = c("Bias", "Coverage probability",
                   "Coverage probability"),
        y = c(0, 0.95, 0.5),
        qntl = c('fake', '95%', '50%')
      )
      dummy$metric <- factor(
        dummy$metric, levels = levels(y$metric)
      )
      dummy2$metric <- factor(
        dummy2$metric, levels = levels(y$metric)
      )
      p <- ggplot(y) +
        geom_point(
          aes(true_eps, med, col = scenario_type, group = qntl),
          position = position_dodge(width = dodge_width),
          size = 1.2
        ) +
        geom_linerange(
          aes(true_eps, ymin = low, ymax = high,
              col = scenario_type, group = qntl),
          position = position_dodge(width = dodge_width)
        ) +
        geom_blank(
          data = dummy, aes(y = low)
        ) +
        geom_blank(
          data = dummy, aes(y = high)
        ) +
        geom_hline(
          data = dummy2, aes(yintercept = y, group = qntl),
          linetype = "dashed"
        ) +
        facet_wrap(~metric, scales = "free_y", ncol = 2) +
        theme_manuscript() +
        scale_color_manual(
          values = values,
          breaks = breaks
        ) +
        xlab("Transmission Advantage") +
        ylab("")
      if (scenario == "same_si") {
        p <- p + theme(legend.position = "none")
      }
      if (scenario == "underrep") {
        p <- p + labs(color = "Reporting probability")
      } else {
        p <- p + labs(color = "Scenario Type")
      }
      save_multiple(
        p, glue("figures/{scenario}_{index}")
      )
    }
    )
  }
)

