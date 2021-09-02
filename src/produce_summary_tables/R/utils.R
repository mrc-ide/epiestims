multiplier_label <- function(val, ref) {
  paste("X", round(val/ref, 1))
}

theme_manuscript <- function(base_size = 14) {
  theme_minimal() %+replace%
    theme(
      text = element_text(size = base_size),
      legend.position = "top"
    )
}
## give filename without the extension
save_multiple <- function(plot, filename) {
  ggsave(
    filename = glue("{filename}.pdf"),
    plot
  )
  ggsave(
    filename = glue("{filename}.png"),
    plot)
}


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

rt_labeller <- function(val) {
  paste("Reference Rt:", val)
}

tmax_labeller <- function(val) {
  paste(val, "days")
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
    theme_manuscript()
  p
}

classification_fig <- function(df) {
  p <- ggplot(df) +
  geom_point(
    aes(as.factor(true_eps), val, col = label)
  ) +
  scale_color_discrete(
    breaks = c("low_greater_than_1",
               "high_less_than_1",
               "CrI_includes_1"),
    labels = c("More transmissible",
               "Less transmissible",
               "Unclear")
  ) +
  facet_wrap(
    ~tmax, labeller = labeller(tmax = tmax_labeller),
    ncol = 2
  ) +
    ylim(0, 1) +
  xlab("True transmission advantage") +
  ylab("Proportion") +
  theme_manuscript() +
  theme(legend.title = element_blank())
  p
}
