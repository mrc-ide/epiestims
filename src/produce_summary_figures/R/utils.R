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

true_epsilon_vs_error <- function(x) {
  p <- ggplot(x) +
    geom_point(
      aes(true_eps, med, col = rt_ref),
      position = position_dodge(width = dodge_width),
      size = 2
    ) +
    geom_linerange(
      aes(true_eps, ymin = low, ymax = high, col = rt_ref),
      position = position_dodge(width = dodge_width),
      size = 1.1
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Estimated - True transmission advantage") +
    xlab("True transmission advantage") +
    labs(color = "Reference Rt") +
    theme_manuscript()
  p
}
