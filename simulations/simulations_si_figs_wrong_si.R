nearest_10 <- function(x, base = 10) 10^ceiling(log(x, base = base))
incid_at_tmax <- readRDS(
  "results/incid_at_tmax_wrong_si.rds"
)

eps_err_summary <- readRDS(
  "results/eps_err_summary_wrong_si.rds"
)

eps_summary <- readRDS(
  "results/eps_summary_wrong_si.rds"
)

## Number of simulations where true epsilon is in
## 95% CrI
by_tmax <- group_by(eps_summary, tmax, epsilon) %>%
  summarise(
    n = sum(epsilon > `2.5%` & epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()

p <- ggplot(by_tmax, aes(tmax, n / total)) +
  geom_point() +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "top")

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
  geom_point(aes(ratio_incid_round, n / total)) +
  ##scale_x_continuous(limits = c(NA, 15)) +
  facet_wrap(~epsilon) +
  ylim(0, 1) +
  theme_minimal()

ggsave("figures/by_tmax_prop_wrong_si.pdf", p)
ggsave("figures/by_ref_incid_wrong_si.pdf", p1)
ggsave("figures/by_var_incid_wrong_si.pdf", p2)
ggsave("figures/by_ratio_wrong_si.pdf", p3)

