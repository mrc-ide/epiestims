nearest_10 <- function(x, base = 10) 10^ceiling(log(x, base = base))
prefix <- "vary_si"
incid_at_tmax <- readRDS(glue("results/incid_at_tmax_{prefix}.rds"))

eps_err_summary <- readRDS(glue("results/eps_err_summary_{prefix}.rds"))

eps_summary <- readRDS(
  glue("results/eps_summary_{prefix}.rds")
)

by_si <- group_by(eps_summary, si_mu_variant, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()
by_si$label <- glue("X {round(by_si$si_mu_variant / 6.83, 1)}")

p <- ggplot(by_si) +
  geom_point(aes(label, n / total)) +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("figures/by_sivar_prop_{prefix}.png", p)

## Number of simulations where true epsilon is in
## 95% CrI
by_tmax <- group_by(eps_summary, tmax, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n(),
    lower = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 2],
    upper = Hmisc::binconf(x = n, n = total, alpha = 0.05)[1, 3]
  ) %>% ungroup()

p <- ggplot(by_tmax) +
  geom_linerange(aes(x = tmax, ymin = lower, ymax = upper)) +
  geom_point(aes(tmax, n / total), size = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(glue("figures/by_tmax_prop_{prefix}.png"), p)
##
eps_summary_incid <- left_join(eps_summary, incid_at_tmax)
## As we will have very different values of incid
## we can summarise usefully by rounding it to the
## nearest power of 10.
eps_summary_incid$ref_incid_round <- nearest_10(eps_summary_incid$ref_incid)
eps_summary_incid$var_incid_round <- nearest_10(eps_summary_incid$var_incid)
eps_summary_incid$ratio_incid_round <- nearest_10(eps_summary_incid$var_incid / eps_summary_incid$ref_incid)
## Some crazy-high numbers here. Restrict to sane
## values
eps_summary_incid <- eps_summary_incid[eps_summary_incid$ratio_incid_round <= 100, ]

by_ref_incid <- group_by(eps_summary_incid, ref_incid_round, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()

by_var_incid <- group_by(eps_summary_incid, var_incid_round, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()

by_ratio <- group_by(eps_summary_incid, ratio_incid_round, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()

p1 <- ggplot(by_ref_incid) +
  geom_point(aes(ref_incid_round, n / total)) +
  scale_x_log10() +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal()


p2 <- ggplot(by_var_incid) +
  geom_point(aes(log(var_incid_round), n / total)) +
  ##scale_x_continuous(limits = c(NA, 15)) +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal()


p3 <- ggplot(by_ratio) +
  geom_point(aes(log(ratio_incid_round), n / total)) +
  ##scale_x_continuous(limits = c(NA, 15)) +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal()


ggsave(glue("figures/by_ref_incid_{prefix}.png"), p1)
ggsave(glue("figures/by_var_incid_{prefix}.png"), p2)
ggsave(glue("figures/by_ratio_{prefix}.png"), p3)

#####
eps_err_summary_incid <- left_join(eps_err_summary, incid_at_tmax)


p <- ggplot(eps_err_summary, aes(factor(tmax), `50%`)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "top")


p <- ggplot(eps_err_summary_incid, aes(log(ref_incid), `50%`)) +
  geom_point() +
  theme_minimal() +
  xlim(0, 15) +
  theme(legend.position = "top")

large <- eps_err_summary_incid[eps_err_summary_incid$`50%` > 10, ]
## cases when error is smalll
small <- eps_err_summary_incid[eps_err_summary_incid$`50%` < 10, ]

p <- ggplot(small, aes(factor(si_mu_variant), `50%`)) +
  geom_boxplot() +
  theme_minimal()

ggsave(glue("figures/by_sivar_{prefix}.png"), p)


p <- ggplot(by_tmax, aes(tmax, epsilon, fill = n)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd")

ggsave(glue("figures/by_tmax_prop_{prefix}_heatmap.png"), p)
