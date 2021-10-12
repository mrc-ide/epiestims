##orderly::orderly_develop_start()
source("R/fig_utils.R")
dir.create("figures")
palette <- c(
  wildtype = "#0f0e0e",
  alpha = "#E69F00",
  betagamma = "#56B4E9",
  delta = "#009E73"
)

variant_nicenames <- c(
  wildtype = "Wildtype", alpha = "Alpha",
  delta = "Delta", betagamma = "Beta/Gamma"
)

## Panel A. Plot of incidence of variants
incidence <- readRDS("cuml_incid_all_variants.rds")
## Repeat Frnech data once to get betagamma as well.
incidence[["french_betagamma"]] <- incidence[["french"]]

tall_incid <- map2(
  incidence,
  list(
    french = c("wildtype", "alpha"),
    uk_alpha_wild = c("wildtype", "alpha"),
    uk_delta_alpha = c("alpha", "delta"),
    french_betagamma = c("wildtype", "betagamma")
  ),
  function(x, y) {
    y <- c("date", y)
    x <- x[, y]
    out <- gather(x, variant, incid, -date)
    out$date <- as.Date(out$date)
    out
  }
)

incid_plots <- map(
  tall_incid, function(x) {
    variants <- intersect(names(palette), unique(x$variant))
    values <- palette[variants]
    ggplot(x) +
      geom_line(
        aes(date, incid, col = variant),
        size = 1.1
      ) +
      scale_color_manual(
        values = values, labels = variant_nicenames
      ) +
      scale_x_date(
        date_breaks = "2 weeks",
        date_labels = "%d-%b-%Y"
      ) +
      ylab("Daily incidence") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )
  }
)
iwalk(
  incid_plots, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_incidence")
    )
  }
)
