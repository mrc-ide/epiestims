## orderly::orderly_develop_start()
## Aesthetics
dir.create("figures")
dodge_width <- 0.5
## common stuff
ms_tmax <- "50"
## SIs of interest
ms_si <- c("X 0.5", "X 1.2", "X 1.5")
si_mu_ref <- 5.4
si_std_ref <- 1.5

vary_si_err <- readRDS(
  "vary_si_err_summary_by_all_vars.rds"
)
eps_vals <- unique(vary_si_err$true_eps)
vary_si_err$true_eps <- factor(
  vary_si_err$true_eps,
  levels = eps_vals, ordered = TRUE
)
vary_si_err$rt_ref <- factor(
  vary_si_err$rt_ref
)

vary_si_err$label <- multiplier_label(
  vary_si_err$si_mu_variant, si_mu_ref
)

## Separate the case when si_mu_variant = si_ref_variant
same_si_mu <- vary_si_err[vary_si_err$label == "X 1", ]
same_si_mu1 <- same_si_mu[same_si_mu$tmax == ms_tmax, ]
## These are the only ones we want to show, in Main and
## supplementary text
vary_si_err <- vary_si_err[vary_si_err$label %in% ms_si, ]
## Numbers for results section
## Median error across Rt_tef and SI_mu at
## various tmax values. To show that error becomes small.
## the bounds represent the variation in median error.
vary_si_err_tab <- select(vary_si_err, -label) %>%
  group_by(true_eps, tmax) %>%
  summarise(
    median_low = quantile(med, 0.025),
    median_med = quantile(med, 0.5),
    median_high = quantile(med, 0.975)
  ) %>% ungroup()
## Format
vary_si_err_tab <- mutate_if(
  vary_si_err_tab, is.numeric, round, 3
)
vary_si_err_tab$formatted <-
  glue("{vary_si_err_tab$median_med}",
       " ({vary_si_err_tab$median_low}, {vary_si_err_tab$median_high})"
       )
vary_si_err_tab <- select(vary_si_err_tab, true_eps, tmax, formatted) %>%
  spread(key = tmax, value = formatted)
## Ordered factors are being converted to
## integers when dumping them into a file.
vary_si_err_tab$true_eps <- eps_vals
cat(
  stargazer(
    vary_si_err_tab[, c("true_eps", "10", "50")], summary = FALSE
  ),
  file = "vary_si_error.tex"
)

## Main text figures
vary_si_err1 <- vary_si_err[vary_si_err$tmax == ms_tmax, ]

p1a <- true_epsilon_vs_error(same_si_mu1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    ) +
  theme(legend.position = "none")
save_multiple(p1a, "figures/same_si_error")

p1b <- true_epsilon_vs_error(vary_si_err1, "Variant SI Mean") +
      facet_wrap(
      ~rt_ref, ncol = 1,
      labeller = labeller(rt_ref = rt_labeller)
    )

save_multiple(p1b, "figures/vary_si_error")

## Supplementary figures
psi <- true_epsilon_vs_error(vary_si_err, "Variant SI Mean") +
      facet_grid(
      tmax~rt_ref,
      labeller = labeller(rt_ref = rt_labeller,
                          tmax = tmax_labeller)
    )
save_multiple(psi, "figures/vary_si_error_by_tmax")


vary_si_classified <- readRDS("vary_si_classified.rds")

p <- ggplot(vary_si_classified) +
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
  xlab("True transmission advantage") +
  ylab("Proportion") +
  theme_manuscript() +
  theme(legend.title = element_blank())

save_multiple(p, "figures/vary_si_classification")
