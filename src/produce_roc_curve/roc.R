## orderly::orderly_develop_start(use_draft = "newer")
source("R/classify_utils.R")
summary_tmax_eps <- function(x) {
  ## This method has the advantage
  ## that numbers across the three
  ## estimated class will add to the
  ## total number of simulations
  out <- group_by(x, across(!sim & !est_class)) %>%
    summarise(
      `Variant more transmissible` = sum(est_class == "Variant more transmissible"),
      `Not more transmissible/Unclear` = sum(est_class == "Not more transmissible/Unclear")
    ) %>% ungroup()
  ## This will of course be the number of sims
  out$nsims <- out$`Not more transmissible/Unclear` +
    out$`Variant more transmissible`
  out
}

infiles <- list(
  vary_offs = "vary_offs_eps_summary_df.rds"
)

eps_summary <- map(infiles, readRDS)
probs <- seq(from = 0.025, to = 0.975, by = 0.025)
low <- scales::label_percent()(probs)

classified <- map(
  eps_summary, function(x) {
    x$true_label <- true_class(x)
    map_dfr(low, function(thres) {
      x$est_class <- "Not more transmissible/Unclear"
      x$est_class[x[[thres]] > 1] <- "Variant more transmissible"
      x <- select(
        x, tmax, sim, rt_ref:est_class
      )
      x$threshold <- thres
      x
    })
  })

x <- summary_tmax_eps(classified[[1]])




x$label <- case_when(
  x$true_label == "Variant more transmissible" &
  x$est_class == "Variant more transmissible" ~ "True Positive",

  x$true_label == "Variant more transmissible" &
  x$est_class == "Not more transmissible/Unclear" ~ "False Negative",

  x$true_label == "No transmission advantage" &
  x$est_class == "Not more transmissible/Unclear" ~ "True Negative",

  x$true_label == "No transmission advantage" &
  x$est_class == "Variant more transmissible" ~ "False Positive"
)

## This makes our x-axis
eps1 <- x[x$true_label == "No transmission advantage", ]
## Can only get specificity from this
eps1$specificity <- eps1$`Not more transmissible/Unclear` / eps1$nsims
## This mkaes our y-axis
## Can only get sensitivity from this
epsy <- x[x$true_label != "No transmission advantage", ]
epsy <- group_by(epsy, tmax, rt_ref, true_eps, kappa, true_label, threshold) |>
  summarise(sensitivity = `Variant more transmissible` / nsims)

## Now put them together
out <- left_join(epsy, eps1, by = c("tmax", "rt_ref", "kappa", "threshold"))

out10 <- out[out$tmax == 50, ]
ggplot(out10, aes(specificity, sensitivity, col = true_eps.x)) + geom_point() +
  facet_grid(rt_ref~kappa)
