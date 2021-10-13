library(dplyr)
library(ggplot2)
eps_qntls <- readRDS("epsilon_qntls_per_region.rds")
## Alpha from France and England
x <- eps_qntls$french
x <- x[x$variant == "alpha_vs_wild", ]
x$country <- "France"
x <- arrange(x, `50%`)
y <- eps_qntls$uk_alpha_wild
y$country <- "England"
y <- arrange(y, `50%`)
z <- rbind(x, y)
z$region <- factor(
  z$region,
  levels = c(x$region, y$region),
  ordered = TRUE
)

p <- ggplot(z) +
  geom_point(aes(region, `50%`, col = country)) +
  geom_linerange(
    aes(region, ymin = `2.5%`, ymax = `97.5%`, col = country)
  ) +
  expand_limits(y = 1) +
  ylab("Effective Transmission Advantage of Alpha") +
  coord_flip() +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )
ggsave("figs/other_figs/regional_estimates_alpha.png", p)


## Same for Delta in England
y <- eps_qntls$uk_delta_alpha
y$country <- "England"
y <- arrange(y, `50%`)
y$region <- factor(y$region, levels = y$region, ordered = TRUE)
p <- ggplot(y) +
  geom_point(aes(region, `50%`, col = country)) +
  geom_linerange(
    aes(region, ymin = `2.5%`, ymax = `97.5%`, col = country)
  ) +
  expand_limits(y = 1) +
  ylab("Effective Transmission Advantage of Delta") +
  coord_flip() +
  scale_color_manual(
    values = c(England = "#0072B2")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )

ggsave("figs/other_figs/regional_estimates_delta.png", p)

## Finally Beta/Gamma
x <- eps_qntls$french
x$country <- "France"
x <- x[x$variant != "alpha_vs_wild", ]
x <- arrange(x, `50%`)
x$region <- factor(x$region, levels = x$region, ordered = TRUE)
p <- ggplot(x) +
  geom_point(aes(region, `50%`, col = country)) +
  geom_linerange(
    aes(region, ymin = `2.5%`, ymax = `97.5%`, col = country)
  ) +
  expand_limits(y = 1) +
  ylab("Effective Transmission Advantage of Beta/Gamma") +
  coord_flip() +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )
ggsave("figs/other_figs/regional_estimates_beta-gamma.png", p)
