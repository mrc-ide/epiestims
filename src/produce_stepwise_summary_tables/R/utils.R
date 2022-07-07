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
    ylab("Proportion classified as") +
    theme_manuscript() +
    theme(legend.title = element_blank())
  p
}
