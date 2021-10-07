library(EpiEstim)
library(glue)
library(purrr)
## estimates by week
## when is the earliest we could have estimated
## epsilon
cbind_rep <- function(x, n) {
  matrix(x, nrow = length(x), ncol = n, byrow = FALSE)
}
mypercent <- function(x) scales::percent(x, accuracy = 0.1)
with_prop_variant <- function(x) {
  p <- ggplot(x) +
    geom_point(aes(col_prop, `50%`)) +
    geom_linerange(aes(x = col_prop, ymin = `2.5%`, ymax = `97.5%`)) +
    expand_limits(y = c(0, NA)) +
    geom_hline(
      yintercept = 1, linetype = "dashed"
    ) +
    scale_x_continuous(
      labels = mypercent
    ) +
    ylab("Effective Transmission Advantage") +
    theme_minimal() +
    theme(text = element_text(size = 16))
  p
}

t_min <- 14L
priors <- default_priors()
mcmc_controls <- list(
  n_iter = 20000L, burnin = 5000L, thin = 10L
)
epi_params <- readRDS('Rdata/Epi_param.rds')

infiles <- list(
  french = 'Rdata/I_fr.rds',
  uk_alpha_wild = 'Rdata/I_UK1.rds',
  uk_delta_alpha = 'Rdata/I_UK2.rds'
)

incidence <- map(infiles, readRDS)

incid_array <- map(
  incidence, function(x) {
    time <- nrow(x[[1]])
    nlocation <- ncol(x[[1]]) - 1 ## -1 for date
    nvariant <- length(x)
    incid <- array(
      NA, dim = c(time, nlocation, nvariant)
    )
    for(i in 1:nvariant) {
      incid[,,i] <- as.matrix(x[[i]][, -1])
    }
  incid
})


estimates <- map2(
  incid_array, incidence, function(x, df) {
    t_max <- seq(
      from = t_min + 7,
      to = dim(x)[1], by = 7
    )
    out <- map(
      t_max, function(tmax) {
        message("t_max = ", tmax)
        estimate_joint(
          incid = x,
          si_distr = cbind_rep(x = Epi_param$SI, n = dim(x)[3]),
          mcmc_control = mcmc_controls,
          priors = priors,
          t_min = t_min,
          t_max = as.integer(tmax)
        )
    })
    names(out) <- df[[1]][["date"]][t_max]
    out
  }
)

eps_estimates <- map2(
  estimates,
  list(
    french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
    uk_alpha_wild = c("alpha_vs_wild"),
    uk_delta_alpha = c("delta_vs_alpha")
  ),
  function(region, variants) {
    map_dfr(region, function(x) {
      out <- apply(
        x[["epsilon"]], 1, quantile,
        prob = c(0.025, 0.5, 0.975)
      )
      out <- data.frame(out)
      names(out) <- variants
      out <- tibble::rownames_to_column(out, "qntl")
      out <- tidyr::gather(out, variant, epsilon, -qntl)
      tidyr::spread(out, qntl, epsilon)
    }, .id = "date")
  })

advantage_with_time <- function(y, var) {
  p <- ggplot(y) +
    geom_point(aes(date, `50%`)) +
    geom_linerange(aes(x = date, ymin = `2.5%`, ymax = `97.5%`)) +
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

volzetal <- data.frame(
  date = "2020-12-31",
  ymin = 1.4, ymax = 1.8
)

p <- ggplot() +
  geom_point(data = alpha_estimates, aes(date, `50%`)) +
  geom_linerange(
    data = alpha_estimates, aes(x = date, ymin = `2.5%`, ymax = `97.5%`)
  ) +
  aes(col = region) +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  geom_linerange(
    data = volzetal, aes(date, ymin = ymin, ymax = ymax),
    col = "black"##, position = position_nudge(x = 0.2)
    ##inherit.aes = FALSE
  ) +
  expand_limits(y = c(0, NA)) +
  geom_hline(
    yintercept = 1, linetype = "dashed"
  ) +
  ylab("Effective Transmission Advantage of Alpha") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    legend.position = "bottom", legend.title = element_blank()
  )
ggsave("figs/other_figs/both_epsilon_alpha.png", p)
## National incidence
afr_total_incid <- data.frame(
  date = incidence[["french"]][["wild"]][["date"]],
  wildtype = apply(incidence[["french"]][["wild"]][, -1], 1, sum),
  alpha = apply(incidence[["french"]][["alpha"]][, -1], 1, sum),
  betagamma = apply(incidence[["french"]][["beta/gamma"]][, -1], 1, sum)
)

uk1_total_incid <- data.frame(
  date = incidence[["uk_alpha_wild"]][["wild"]][["date"]],
  wildtype = apply(incidence[["uk_alpha_wild"]][["wild"]][, -1], 1, sum),
  alpha = apply(incidence[["uk_alpha_wild"]][["alpha"]][, -1], 1, sum)
)

uk2_total_incid <- data.frame(
  date = incidence[["uk_delta_alpha"]][["alpha"]][["date"]],
  alpha = apply(incidence[["uk_delta_alpha"]][["alpha"]][, -1], 1, sum),
  delta = apply(incidence[["uk_delta_alpha"]][["delta"]][, -1], 1, sum)
)



cuml_incid <- map(
  list(
    french = fr_total_incid,
    uk_alpha_wild = uk1_total_incid,
    uk_delta_alpha = uk2_total_incid
  ), function(x) {
    out <- apply(x[, -1], 2, cumsum)
    prop_variant <- out / apply(out, 1, sum)
    res <- cbind(x, out, prop_variant)
    names(res) <- c(
      names(x), glue("cumulative_{names(x[, -1])}"),
      glue("proportion_{names(x[ ,-1])}")
    )
    res
  }
)

## x-axis is now the proportion of the variant cases
joined <- map2(
  eps_estimates, cuml_incid, function(x, y) {
    x$date <- as.Date(x$date)
    left_join(x, y, by = "date")
  }
)
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
  xlab("Proportion of Alpha") +
  scale_color_manual(
    values = c(England = "#0072B2", France = "#D55E00")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
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

