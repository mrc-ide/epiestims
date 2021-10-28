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
    xlab("True Transmission Advantage") +
    ylim(0, 1) +
    theme_manuscript() +
    labs(color = "Reference Rt") +
    theme(legend.position = "top")
  p
}


true_epsilon_vs_error <- function(x, color_by) {
  p <- ggplot(x) +
    geom_point(
      aes(true_eps, med), col = "black",
      position = position_dodge(width = dodge_width),
      size = 2
    ) +
    geom_linerange(
      aes(true_eps, ymin = low, ymax = high), col = "black",
      position = position_dodge(width = dodge_width),
      size = 1.1
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Bias") +
    xlab("True Transmission Advantage") +
    labs(color = color_by) +
    theme_manuscript(base_size = 14)
  p
}

classification_fig <- function(df) {
  p <- ggplot(df) +
    geom_point(
      aes(true_eps, PointEst), col = "black",
      size = 1.4
    ) +
    geom_linerange(
      aes(true_eps, ymin = `Lower`, ymax = `Upper`), col = "black",
      size = 1.1
    ) +
    facet_grid(
      tmax~rt_change,
      labeller = labeller(tmax = tmax_labeller,
                          rt_change = rt_change_labeller)
    ) +
    xlab("True Transmission Advantage") +
    ylab("Proportion classified correctly") +
    theme_manuscript(base_size = 14) +
    theme(legend.title = element_blank())
  p
}
