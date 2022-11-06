## orderly::orderly_develop_start(use_draft = "newer")
source("R/classify_utils.R")
source("R/fig_utils.R")
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
  baseline = "vary_si_eps_summary_df.rds",
  vary_offs = "vary_offs_eps_summary_df.rds"
)

eps_summary <- map(infiles, readRDS)
si_mu_ref <- 5.4
si_std_ref <- 1.5
x <- eps_summary[["baseline"]]
x <- x[x$si_mu_variant == si_mu_ref, ]
eps_summary[["baseline"]] <- x

low <- grep(pattern = "%", x = colnames(eps_summary[[1]]), value = TRUE)

classified <- map(
  eps_summary, function(x) {
    x$true_label <- true_class(x)
    map_dfr(low, function(thres) {
      x$est_class <- "Not more transmissible/Unclear"
      x$est_class[x[[thres]] > 1] <- "Variant more transmissible"
      x <- select(x, tmax, sim, rt_ref:est_class)
      x$threshold <- thres
      summary_tmax_eps(x)
    })
  })

classified[["baseline"]]$superspreading <- "Baseline"
classified[["vary_offs"]]$superspreading <- case_when(
  classified[["vary_offs"]]$kappa == 0.1 ~ "High",
  classified[["vary_offs"]]$kappa == 0.5 ~ "Moderate",
  classified[["vary_offs"]]$kappa == 1 ~ "Low"
)
common <- intersect(colnames(classified[[1]]), colnames(classified[[2]]))
x <- rbind(classified[[1]][, common], classified[[2]][, common])
## This makes our x-axis
eps1 <- x[x$true_label == "No transmission advantage", ]
## Can only get specificity from this
eps1$specificity <- eps1$`Not more transmissible/Unclear` / eps1$nsims
## This mkaes our y-axis
## Can only get sensitivity from this
epsy <- x[x$true_label != "No transmission advantage", ]
epsy <- group_by(epsy, tmax, rt_ref, true_eps, superspreading, true_label, threshold) |>
summarise(sensitivity = `Variant more transmissible` / nsims) |>
ungroup()

## Now put them together
## keep only relevant columns to minimise confusiom
cols <- c("tmax", "rt_ref", "true_eps", "superspreading", "threshold", "sensitivity")
epsy <- epsy[, cols]

cols <- c("tmax", "rt_ref", "superspreading", "threshold", "specificity")
eps1 <- eps1[, cols]

out <- left_join(epsy, eps1)
out$fpr <- 1 - out$specificity
out <- out[out$`true_eps` %in% c(1.1, 1.5, 3), ]
out <- split(out, out$rt_ref)

## x <- out[[1]]
## x$threshold_numeric <- readr::parse_number(x$threshold)
## x <- pivot_longer(
##   x, cols = c("sensitivity", "fpr"), names_to = "var"
## )
## ##x <- x[x$tmax == 10, ]

## ggplot(x) +
##   geom_point(aes(threshold_numeric, fpr, col = true_eps)) +
##   facet_grid(tmax~superspreading) +
##   theme_manuscript() +
##     theme(legend.position = "top")


iwalk(out, function(x, rt) {

  x$superspreading <- factor(
    x$superspreading, levels = c("Baseline", "Low", "Moderate", "High"),
    ordered = TRUE
  )
  y <- x[x$threshold == "2.5%", ]
  p <- ggplot() +
    geom_line(data = x, aes(fpr, sensitivity, col = `tmax`)) +
    geom_point(
      data = y, aes(fpr,sensitivity, col = tmax), shape = 8,
      size = 2
    ) +
  facet_grid(`true_eps`~superspreading) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    scale_colour_viridis_d(
      option = "viridis", name = "Days of data used"
    ) +
      theme_manuscript() +
    theme(legend.position = "top")

  save_multiple(p, glue("roc_{rt}"))
})

