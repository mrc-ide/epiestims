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
  
  together <- joined_data
  x1 <- together[[1]]
  x1$qntl <- 'fake'
  x2 <- together[[2]]
  x2$qntl <- 'fake'
  x3 <- together[[3]]
  ## This now has 50% coverage probability as well
  x31 <- select(x3, rt_ref:upper, rt_change:metric)
  x32 <- select(
    x3, rt_ref:true_eps, n = n50, total = total, pt_est = pt_est50,
    lower = lower50, upper = upper50, rt_change:metric
  )
  x31$qntl <- '95%'
  x32$qntl <- '50%'
  x3 <- rbind(x31, x32)
  x3 <- rename(
    x3, low = lower, med = pt_est, high = upper
  )
  x3$si_mu_variant <- 5.4
  x3$label <- "X 1"
  x3 <- x3[, colnames(x2)]
  x4 <- together[[4]]
  x4$qntl <- '1'
  idx1 <- which(x4$true_label == "No transmission advantage" &
                  x4$est_class == "Unclear")
  idx2 <- which(x4$true_label == x4$est_class)
  # browser()
  x4 <- x4[c(idx1, idx2), ]
  x4$label <- "X 1"
  # x4 <- rename(
  #   x4, med = PointEst, low = Lower, high = Upper
  # )
  x4 <- x4[, colnames(x2)]
  out <- rbind(x1, x2, x3, x4)
  
  out$metric <- factor(
    out$metric, levels = c("Bias", "Uncertainty",
                                "Coverage probability",
                                "Classification"),
    ordered = TRUE
  )
  
  ## Construct dummy data.frame to control facet scales
  ## Coverage probability to go from 0 to 1
  min_bias <- -1.5
  max_bias <- 1.5
  min_sd <- -0.1
  max_sd <- 0.75
  dummy <- data.frame(
    metric = c("Bias", "Coverage probability", "Uncertainty", "Classification"),
    ##true_eps = levels(y$true_eps),
    low = c(min_bias, 0, min_sd, 0),
    high = c(max_bias, 1, max_sd, 1)
  )
  dummy2 <- data.frame(
    metric = c("Bias", "Coverage probability",
               "Coverage probability"),
    y = c(0, 0.95, 0.5),
    qntl = c('fake', '95%', '50%')
  )
  
  dummy$metric <- factor(
    dummy$metric, levels = levels(out$metric)
  )
  dummy2$metric <- factor(
    dummy2$metric, levels = levels(out$metric)
  )
  
  ## Different shapes for 95% and 50% coverage probability
  out$shape <- 19 ## everything is a circle
  out$shape[out$qntl == "50%"] <- 18 ## Except 50% coverage probability
  
  
  y <- split(out, list(out$tmax, out$rt_change))
  iwalk(y, function(z, index) {
    
    p <- ggplot(z) +
      geom_point(
        aes(true_eps, med, group = qntl, shape = shape), col = "black",
        position = position_dodge2(width = dodge_width),
        size = 1.2
      ) +
      geom_linerange(
        aes(true_eps, ymin = low, ymax = high, group = qntl),
        col = "black",
        position = position_dodge2(width = dodge_width)
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
      scale_shape_identity() +
      theme_manuscript() +
      xlab("Transmission Advantage") +
      ylab("")
    
    save_multiple(
      p, glue("figures/{panel_name}_panel_{index}")
    )
    
  })
  
  
  
}
