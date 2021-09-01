## orderly::orderly_develop_start()
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

## Aesthetics
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"

si_mu_ref <- 5.4
si_std_ref <- 1.5

vary_si_err <- readRDS(
  "vary_si_err_summary_by_all_vars.rds"
)

vary_si_err$si_label <- multiplier_label(
  vary_si_err$si_mu_variant, si_mu_ref
)

eps_vals <- unique(vary_si_err$true_eps)
vary_si_err$true_eps <- factor(
  vary_si_err$true_eps,
  levels = eps_vals, ordered = TRUE
)
vary_si_err$rt_ref <- factor(
  vary_si_err$rt_ref
)

## Separate the case when si_mu_variant = si_ref_variant
same_si_mu <- vary_si_err[vary_si_err$si_label == "X 1", ]
vary_si_err <- vary_si_err[vary_si_err$si_label != "X 1", ]


## Main text figures
same_si_mu1 <- same_si_mu[same_si_mu$tmax == ms_tmax, ]
vary_si_err1 <- vary_si_err[vary_si_err$tmax == ms_tmax, ]

p1a <- true_epsilon_vs_error(same_si_mu1)
save_multiple(p1a, "figures/fig1a")


p1b <- true_epsilon_vs_error(vary_si_err1) +
  facet_wrap(~si_label)
