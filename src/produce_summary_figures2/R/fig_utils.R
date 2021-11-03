scenario_type_labeller <- function(x) {
  x$scenario_type <- case_when(
    x$label == "X 0.5" ~ "Low",
    x$label == "X 1.5" ~ "Moderate",
    x$label == "X 1" ~ "Baseline",
    x$label == "X 2" ~ "High",
    ## Vary offspring
    x$label == "0.1" ~ "High",
    x$label == "0.5" ~ "Moderate",
    x$label == "1" ~ "Low",
    ## Underreporting
    x$label == "0.2" ~ "High",
    x$label == "0.5" ~ "Moderate",
    x$label == "0.8" ~ "Low"
  )
  x
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


multiplier_label <- function(val, ref) {
  paste("X", round(val/ref, 1))
}

theme_manuscript <- function(base_size = 20) {
  theme_classic() %+replace%
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
