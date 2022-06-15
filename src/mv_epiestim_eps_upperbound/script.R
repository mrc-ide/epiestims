## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
eps_upperbound <- function(x) 1/sqrt(x)

incid <- 1:1000
cv <- eps_upperbound(incid)

p <- ggplot() +
  geom_line(aes(incid, cv)) +
  theme_manuscript() +
  xlab("Cumulative Incidence") +
  ylab("Effective Transmission Adavntage CV")

save_multiple(p, 'eps_upperbound')


## Realised CV, not a useful comparison
## nonovl_est <- readRDS("nonoverlapping_epsilon_estimates.rds")

## nonovl_est_cv <- map2(
##   nonovl_est,
##   list(
##     french = c("alpha_vs_wild", "beta-gamma_vs_wild"),
##     uk_alpha_wild = c("alpha_vs_wild"),
##     uk_delta_alpha = c("delta_vs_alpha")
##   ),
##   function(est, names) {
##     map_dfr(
##       est, function(x) {
##         eps <- x[["epsilon"]]
##         mu <- apply(eps, 1, mean)
##         sigma <- apply(eps, 1, sd)
##         cv <- sigma / mu
##         out <- data.frame(cv = cv)
##         out$var <- names
##         out
##       }, .id = "date"
##     )
##   }
## )

## x <- nonovl_est_cv[["french"]]
## nonovl_est_cv[["french_betagamma"]] <- x[x$var == "beta-gamma_vs_wild", ]
## nonovl_est_cv[["french"]] <- x[x$var != "beta-gamma_vs_wild", ]

## nonovl_prop <- readRDS("nonoverlapping_prop_variant.rds")
## nonovl_prop[["french_betagamma"]] <- nonovl_prop[["french"]]
## nonovl_prop <- map_depth(nonovl_prop, 2, slice_tail)
## nonovl_prop <- map(nonovl_prop, bind_rows)

## ## Theoretical upperbound vs realised CV
## theort_actual <- pmap(
##   list(
##     x = nonovl_est_cv,
##     y = nonovl_prop,
##     cols = list(
##       french = c("cumulative_wildtype", "cumulative_alpha"),
##       uk_alpha_wild = c("cumulative_wildtype", "cumulative_alpha"),
##       uk_delta_alpha = c("cumulative_alpha", "cumulative_delta"),
##       french_betagamma = c("cumulative_wildtype", "cumulative_betagamma")
##     )), function(x, y, cols) {
##       y$total <- y[[cols[1]]] + y[[cols[2]]]
##       y$upper_bound <- eps_upperbound(y$total)
##       x$date <- as.Date(x$date)
##       y$date <- as.Date(y$date)
##       left_join(x, y, by = "date")
##     })

