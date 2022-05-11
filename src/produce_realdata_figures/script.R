## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")
mypercent <- function(x) scales::percent(x)
palette <- c(
  wildtype = "#0f0e0e",
  alpha = "#E69F00",
  betagamma = "#CC79A7",
  `beta/gamma` = "#CC79A7",
  delta = "#009E73",
  England = "#cc0000",
  France = "#0000ff"
)

date_breaks <- "6 weeks"
date_labels <- "%d-%b"


variant_nicenames <- c(
  wildtype = "Wildtype", alpha = "Alpha",
  delta = "Delta", betagamma = "Beta/Gamma",
  `beta/gamma` = "Beta/Gamma"
)

## Read all data first to get consistent ylimits
incidence <- readRDS("cuml_incid_all_variants.rds")
## Repeat Frnech data once to get betagamma as well.
incidence[["french_betagamma"]] <- incidence[["french"]]
date_min <- map(incidence, ~ min(.$date))

xaxis_breaks <- list(
  french = as.Date(c("2021-02-15", "2021-03-15", "2021-04-15", "2021-05-15")),
  uk_alpha_wild = as.Date(
    c("2020-09-01", "2020-10-01", "2020-11-01", "2020-12-01", "2021-01-01",
      "2021-02-01", "2021-03-01"
      )),
  uk_delta_alpha = as.Date(
    c("2021-03-15", "2021-04-15", "2021-05-15", "2021-06-15")
  ),
  french_betagamma = as.Date(
    c("2021-02-15", "2021-03-15", "2021-04-15", "2021-05-15")
  )
)

## Alpha multiplicative advantage by Volz et al
volzetal <- data.frame(
  date = as.Date("2020-12-31"),
  ymin = 1.4, ymax = 1.8, y = 1.74
)

##eps_over_time <- readRDS("epsilon_qntls_over_time.rds")
eps_over_time <- readRDS("nonoverlapping_epsilon_qntls.rds")
eps_over_time[["french_betagamma"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant != "alpha_vs_wild", ]
eps_over_time[["french"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant == "alpha_vs_wild", ]


eps_over_time_with_prop <- readRDS("epsilon_estimates_with_variant_proportion.rds")
eps_over_time_with_prop[["french_betagamma"]] <-
  eps_over_time_with_prop[["french"]][eps_over_time_with_prop[["french"]]$variant != "alpha_vs_wild", ]
eps_over_time_with_prop[["french"]] <-
  eps_over_time_with_prop[["french"]][eps_over_time_with_prop[["french"]]$variant == "alpha_vs_wild", ]

both_together <- pmap(
  list(
    default = eps_over_time_with_prop,
    ##custom = custom_eps_over_time,
    column = c(
      "proportion_alpha", "proportion_alpha", "proportion_delta",
      "proportion_betagamma")
  ), function(default, custom, column) {
    ##default$prior <- "Default prior"
    ##custom$prior <- "Informative prior"
    ##x <- rbind(default[, colnames(custom)], custom)
    x <- default
    x$proportion <- x[[column]]
    x
  }
)

######################################################################
## Plot with 2 y-axis
plots2axis <- pmap(
  list(
    x = eps_over_time,
    y = both_together,
    y2label =  c("Cumulative proportion of Alpha",
                 "Cumulative proportion of Alpha",
                 "Cumulative proportion of Delta",
                 "Cumulative proportion of Beta/Gamma"),
    xmin = xaxis_breaks,
    col = palette[c("alpha", "alpha", "delta", "betagamma")]
  ),
  function(x, y, y2label, xmin, col) {
    y <- select(y, date, proportion)
    x$date <- as.Date(x$date)
    z <- left_join(x, y, by = "date")
    ## coeff <-  max(z$`97.5%`) / max(z$proportion)
    coeff <- 1.8 ## For everyhing except delta
    if (x$variant[1] == "delta_vs_alpha") coeff <- 2.5

    z$proportion_scaled <- z$proportion * coeff
    message("Max z$`97.5%`", max(z$`97.5%`))
    message("Range of proportion", range(z$proportion))
    message("Coeff = ", coeff)
    ggplot(z, aes(x = date)) +
      geom_point(
        aes(y = `50%`), colour = col, size = 3
      ) +
      geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`),
        size = 1.1, colour = col
      ) +
      geom_hline(
        yintercept = 1, linetype = "dashed", color = "red", size = 1.2
      ) +
      geom_line(aes(y = proportion_scaled), color = "blue") +
      scale_y_continuous(
        sec.axis = sec_axis(~./coeff, name = y2label, labels = mypercent),
        limits = c(0, coeff),
        breaks = seq(0, coeff, by = 0.5)
      ) +
      scale_x_date(
        breaks = xmin,
        ## date_breaks = date_breaks,
        date_labels = date_labels,
        limits = c(min(xmin), NA)
      ) +
      coord_cartesian(clip = "off") +
      ylab("Effective transmission advantage") +
      xlab("") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        legend.title = element_blank(),
        axis.line.y.right = element_line(color = "blue"),
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue")
      )

  }
)

iwalk(
  plots2axis, function(p, name) {
    if (name == "uk_alpha_wild") {
      p <- p +
        geom_point(
          data = volzetal, aes(date, y),
          size = 4,
          position = position_nudge(x = 2)
        ) +
        geom_linerange(
          data = volzetal,
          aes(date, ymin = ymin, ymax = ymax),
          size = 1.1,
          position = position_nudge(x = 2)
      )
    }
    save_multiple(p, glue("figures/{name}_2axis"))
  }
)



