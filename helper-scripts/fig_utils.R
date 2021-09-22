multiplier_label <- function(val, ref) {
  paste("X", round(val/ref, 1))
}

theme_manuscript <- function(base_size = 14) {
  theme_minimal() %+replace%
    theme(
      text = element_text(size = base_size),
      legend.position = "top",
      axis.text.x = element_text(angle = 45)
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

rt_labeller <- function(val) {
  paste("Reference Rt:", val)
}

rt_change_labeller <- function(val) {
  paste("Rt change:", val)
}

tmax_labeller <- function(val) {
  paste(val, "days")
}

## df is the output of summarise_median_err
format_median_err <- function(df) {
  df$formatted <-
  glue("{df$median_med}",
       " ({df$median_low}, {df$median_high})")
  df <- select(df, true_eps, tmax, formatted) %>%
    spread(key = tmax, value = formatted)
  df
}

pretty_ci <- function(val, low, high, round_to = 2) {
  f <- function(x) {
    format(round(x, round_to), nsmall = 2)
  }
  glue("{f(val)} \n ({f(low)}, {f(high)})")
}
