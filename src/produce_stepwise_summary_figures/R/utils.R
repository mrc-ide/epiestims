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


true_epsilon_vs_sd <- function(x) {
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
    ylab("Uncertainty") +
    xlab("True Transmission Advantage") +
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

panel_fig <- function(joined_data, panel_name) {
  
  dummy <- data.frame(
    metric = c("Bias", "Coverage probability"),
    ##true_eps = levels(y$true_eps),
    low = 0,
    high = 1
  )
  dummy2 <- data.frame(
    metric = c("Bias", "Coverage probability"),
    y = c(0, 0.95)
  )
  dummy$metric <- factor(
    dummy$metric, levels = levels(joined_data$metric)
  )
  dummy2$metric <- factor(
    dummy2$metric, levels = levels(joined_data$metric)
  )
  
  
  
  y <- split(joined_data, list(joined_data$tmax, joined_data$rt_change))
  iwalk(y, function(z, index) {
    
    p <- ggplot(z) +
      geom_point(
        aes(true_eps, med), col = "black",
        position = position_dodge(width = dodge_width),
        size = 1.2
      ) +
      geom_linerange(
        aes(true_eps, ymin = low, ymax = high),
        col = "black",
        position = position_dodge(width = dodge_width)
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
      xlab("Transmission Advantage") +
      ylab("")
    
    save_multiple(
      p, glue("figures/{panel_name}_panel_{index}")
    )
    
  })
  
  
  
}
