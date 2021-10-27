true_epsilon_vs_95CrI <- function(x) {
  p <- ggplot(x) +
    geom_point(
      aes(true_eps, pt_est, col = rt_ref),
      position = position_dodge(width = dodge_width)
    ) +
    geom_linerange(
      aes(true_eps, ymin = lower, ymax = upper, col = rt_ref),
      position = position_dodge(width = 0.3)
    ) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylab("Proportion in 95% CrI") +
    xlab("True transmission advantage") +
    ylim(0, 1) +
    theme_manuscript() +
    labs(color = "Reference Rt") +
    theme(legend.position = "top")
  p
}


true_epsilon_vs_error <- function(x, color_by) {
  p <- ggplot(x) +
    geom_point(
      aes(true_eps, med, col = label),
      position = position_dodge(width = dodge_width),
      size = 2
    ) +
    geom_linerange(
      aes(true_eps, ymin = low, ymax = high, col = label),
      position = position_dodge(width = dodge_width),
      size = 1.1
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Estimated - True transmission advantage") +
    xlab("True transmission advantage") +
    labs(color = color_by) +
    ##theme_manuscript() +
    theme(legend.position = c(0.7, 0.2))
  p
}

classification_fig <- function(df) {
  p <- ggplot(df) +
  geom_line(
    aes(true_eps, val, col = classification),
    size = 1.2
  ) +
  facet_wrap(
    ~tmax, labeller = labeller(tmax = tmax_labeller),
    ncol = 2
  ) +
  xlab("True transmission advantage") +
  ylab("Proportion") +
  theme_manuscript() +
  theme(legend.title = element_blank())
  p
}

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
 x[["underrep"]]$label <- factor(x[["underrep"]]$p_report)
 x
}

## x is the output of affix_label
main_and_suppl <- function(x, ms_vars, ms_tmax) {
 list(
  same_si = list(
    main = filter(
      x[["vary_si"]], label %in% ms_vars[["same_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["vary_si"]], label %in% ms_vars[["same_si"]],
      tmax != ms_tmax
    )
  ),
  vary_offs = list(
    main = filter(
      x[["vary_offs"]], label %in% ms_vars[["vary_offs"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["vary_offs"]], tmax != ms_tmax
    )
  ),
  vary_si = list(
    main = filter(
      x[["vary_si"]], label %in% ms_vars[["vary_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["vary_si"]], label %in% ms_vars[["vary_si"]],
      tmax != ms_tmax
    )
  ),
  wrong_si = list(
    main = filter(
      x[["wrong_si"]], label %in% ms_vars[["wrong_si"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["wrong_si"]], label %in% ms_vars[["wrong_si"]],
      tmax != ms_tmax
    )
  ),
  vary_cv = list(
    main = filter(
      x[["vary_cv"]], label %in% ms_vars[["vary_cv"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["vary_cv"]], label %in% ms_vars[["vary_cv"]],
      tmax != ms_tmax
    )
  ),
  wrong_cv = list(
    main = filter(
      x[["wrong_cv"]], label %in% ms_vars[["wrong_cv"]],
      tmax == ms_tmax
    ),
    suppl = filter(
      x[["wrong_cv"]], label %in% ms_vars[["wrong_cv"]],
      tmax != ms_tmax
    )
  ),
  underrep = list(
    main = NA, ## All of this goes in the supplementary
    suppl = filter(
      x[["underrep"]], label %in% ms_vars[["underrep"]]
    )
  )
 )
}
## faceting by rt_ref and tmax
suppl_figure <- function(y, index) {
  limits <- intersect(y$scenario_type, names(values))
  y$true_eps <- factor(
    y$true_eps,
    levels = unique(y$true_eps)
  )
   y$scenario_type <- factor(
    y$scenario_type,
    levels = c("Baseline", "Low", "Moderate", "High"),
    ordered = TRUE
  )
  p <- ggplot(y) +
  geom_point(
    aes(true_eps, med, col = scenario_type),
      position = position_dodge(width = dodge_width),
      size = 1.4
  ) +
  geom_linerange(
    aes(true_eps, ymin = low, ymax = high, col = scenario_type),
      position = position_dodge(width = dodge_width),
      size = 1
  ) +
  facet_grid(
    tmax~rt_ref,
    labeller = labeller(tmax = tmax_labeller,
                        rt_ref = rt_labeller)
  ) +
    scale_color_manual(
      values = values,
      breaks = c("Low", "Moderate", "High")
    ) +
  xlab("True Transmssion Advantage") +
  theme_manuscript() +
    theme(legend.position = "bottom")
  p
}
