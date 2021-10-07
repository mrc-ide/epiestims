library(ggplot2)
library(glue)
library(purrr)

mypercent <- function(x) scales::percent(x, accuracy = 0.1)

with_prop_variant <- function(x) {
  p <- ggplot(x) +
    geom_point(aes(col_prop, `50%`)) +
    geom_linerange(aes(x = col_prop, ymin = `2.5%`, ymax = `97.5%`), size = 1.1) +
    expand_limits(y = 1) +
    geom_hline(
      yintercept = 1, linetype = "dashed"
    ) +
    scale_x_continuous(labels = mypercent) +
    ylab("Effective Transmission Advantage") +
    theme_minimal() +
    theme(text = element_text(size = 16))
  p
}

advantage_with_time <- function(y, var) {
  p <- ggplot(y) +
    geom_point(aes(date, `50%`)) +
    geom_linerange(aes(x = date, ymin = `2.5%`, ymax = `97.5%`), size = 1.1) +
    expand_limits(y = c(0, NA)) +
    geom_hline(
      yintercept = 1, linetype = "dashed"
    ) +
    ylab(paste("Effective Transmission Advantage of", var)) +
    theme_minimal() +
    theme(
      text = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)
    )
  p
}

eps_estimates <- readRDS("epsilon_estimates_over_time.rds")
cuml_incid <- readRDS("cuml_incid_all_variants.rds")
joined <- readRDS("epsilon_estimates_with_variant_proportion.rds")

iwalk(
  eps_estimates, function(x, region) {
    variant <- split(x, x$variant)
    iwalk(variant, function(y, var) {
      p <- advantage_with_time(y, var)
      ggsave(
        filename = glue("figs/other_figs/{region}_{var}.png"), p
      )
    })
  }
)
## Can put Alpha estimates from France and England on the same
## figure
eps_estimates[[1]]$region <- "France"
eps_estimates[[2]]$region <- "England"
alpha_estimates <- rbind(
  eps_estimates[[1]][eps_estimates[[1]]$variant == "alpha_vs_wild", ],
  eps_estimates[[2]]
)

volzetal <- data.frame(date = "2020-12-31", ymin = 1.4, ymax = 1.8)

p <- ggplot() +
  geom_point(data = alpha_estimates, aes(date, `50%`)) +
  geom_linerange(
    data = alpha_estimates, aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
    size = 1.1
  ) +
  aes(col = region) +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  geom_linerange(
    data = volzetal, aes(date, ymin = ymin, ymax = ymax), size = 1.1,
    col = "black"##, position = position_nudge(x = 0.2)
    ##inherit.aes = FALSE
  ) +
  expand_limits(y = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylab("Effective Transmission Advantage of Alpha") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    legend.position = "bottom", legend.title = element_blank()
  )
ggsave("figs/other_figs/both_epsilon_alpha.png", p)

## Showing Alpha for England and France on the same graph
x <- joined[["french"]][joined[["french"]]$variant == "alpha_vs_wild", ]
y <- joined[[2]]
x$region <- "France"
y$region <- "England"
## France has more columns than England
z <- rbind(x[, colnames(y)], y)
z$col_prop <- z$proportion_alpha

## On 31st December 2020, when Volz et al was
## published, the proportion of Alpha cases in
## England was:
volzetal <- data.frame(
  proportion_alpha = cuml_incid[[2]][cuml_incid[[2]]$date == "2020-12-31", "proportion_alpha"],
  ymin = 1.4, ymax = 1.8
)


p <- with_prop_variant(z) +
  aes(col = region) +
  geom_linerange(
    data = volzetal,
    aes(proportion_alpha, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE
  ) +
  ylab("Effective Transmission Advantage of Alpha") +
  xlab("Proportion of Alpha") +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.position = "bottom", legend.title = element_blank()
  )

ggsave("figs/other_figs/both_epsilon_proportion_alpha.png", p)
## French Alpha

x <- joined[["french"]][joined[["french"]]$variant == "alpha_vs_wild", ]
x$col_prop <- x$proportion_alpha
p <- with_prop_variant(x) +
  xlab("Proportion of Alpha")
ggsave("figs/other_figs/french_epsilon_proportion_alpha.png", p)

x <- joined[["french"]][joined[["french"]]$variant != "alpha_vs_wild", ]
x$col_prop <- x$proportion_betagamma
p <- with_prop_variant(x) +
  xlab("Proportion of Beta/Gamma")
ggsave("figs/other_figs/french_epsilon_proportion_betagamma.png", p)

x <- joined[[2]]
x$col_prop <- x$proportion_alpha
p <- with_prop_variant(x) +
  xlab("Proportion of Alpha")
ggsave("figs/other_figs/england_epsilon_proportion_alpha.png", p)

x <- joined[[3]]
x$col_prop <- x$proportion_delta
p <- with_prop_variant(x) +
  xlab("Proportion of Delta")
ggsave("figs/other_figs/england_epsilon_proportion_delta.png", p)
