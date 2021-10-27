multiplier_label <- function(val, ref) {
  paste("X", round(val/ref, 1))
}

theme_manuscript <- function(base_size = 20) {
  theme_minimal() %+replace%
    theme(
      text = element_text(size = base_size),
      axis.line = element_line(size = 1.05),
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
