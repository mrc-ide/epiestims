## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")

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




## Panel A. Plot of incidence of variants
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

incid_plots <- map2(
  tall_incid, xaxis_breaks, function(x, breaks) {
    variants <- intersect(names(palette), unique(x$variant))
    values <- palette[variants]
    coeff <- max(x$incid)
    ggplot(x) +
      geom_line(
        aes(date, incid, col = variant),
        size = 1.1
      ) +
      scale_color_manual(
        values = values, labels = variant_nicenames
      ) +
      scale_x_date(
        breaks = breaks,
        date_labels = date_labels
      ) +
      scale_y_continuous(
        labels = scales::label_number(suffix = "K", scale = 1e-3),
        ## This is a dummy secondary axis. Keeping everything else
        ## i.e. axis labels and tick labels the same to help with alignment.
        sec.axis = sec_axis(~./coeff, name = "Cumulative proportion of Alpha",
                            labels = mypercent)

      ) +
      expand_limits(x = range(breaks)) +
      coord_cartesian(clip = "off") +
      ylab("Daily incidence") +
      xlab("") +
  theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        ##Axis.title.x = element_blank(),
        legend.title = element_blank(),
        ##legend.position = c(0.2, 0.85),
        ## We don't actually want to show the right y-axis
        axis.line.y.right = element_line(color = "white"),
        axis.title.y.right = element_text(color = "white"),
        axis.text.y.right = element_text(color = "white"),
        axis.ticks.y.right = element_line(color = "white")
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
## Regional incidence plots for Supplementary
regional_incid <- list(
  french = "I_fr.rds",
  uk1 = "I_UK1.rds",
  uk2 = "I_UK2.rds"
)
regional_incid <- map(regional_incid, readRDS)
regional_incid <- map(
  regional_incid, function(x) bind_rows(x, .id = "variant")
)
regional_incid[["england"]] <- regional_incid[["uk1"]]


walk(
  c("french", "england"),
  function(region) {
    x <- regional_incid[[region]]
    x$variant <- case_when(
      x$variant == "wild" ~ "wildtype",
      TRUE ~ x$variant
    )
    variants <- intersect(names(palette), unique(x$variant))
    values <- palette[variants]

    xtall <- gather(x, region, incid, -date, -variant)
    xtall$date <- as.Date(xtall$date)
    ##xtall$region <- stringr::str_wrap(xtall$region, 25)
    breaks <- seq(
      round_date(min(xtall$date), "month"),
      round_date(max(xtall$date), "month"),
      "1 month"
    )

    p <- ggplot(xtall) +
      geom_line(aes(date, incid, col = variant)) +
      facet_wrap(
        ~region, scales = "free_y", ncol = 3,
        labeller = labeller(region = region_short_names)
      ) +
      scale_color_manual(
        values = values, labels = variant_nicenames,
      ) +
      scale_x_date(
        ##date_breaks = date_breaks,
        breaks = breaks,
        date_labels = date_labels
      ) +
      scale_y_continuous(
        n.breaks = 3
      ) +
      coord_cartesian(clip = "off") +
      ylab("Daily incidence") +
      xlab("") +
      theme_manuscript() +
      theme(
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_blank(),
        legend.position = "top",
        ### Transparent right axis so panels a and c can be aligned
        axis.line.y.right = element_line(color = "blue"),
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue")
      )

    save_multiple(
      p, glue("figures/{region}_regional_incid")
    )
  }
)

########## Plots on log scale ########
tall_incid <- map2(
  incidence,
  list(
    french = c("wildtype", "alpha"),
    uk_alpha_wild = c("wildtype", "alpha"),
    uk_delta_alpha = c("alpha", "delta"),
    french_betagamma = c("wildtype", "betagamma")
  ),
  function(x, y) {
    ydaily <- c("date", y)
    ycum <- c("date", paste0("cumulative_", y))

    x1 <- x[, ydaily]
    ##x1 <- mutate_if(x1, is.numeric, log)
    out <- gather(x1, variant, `Daily Incidence`, -date)
    out$date <- as.Date(out$date)

    x2 <- x[, ycum]
    ##x2 <- mutate_if(x2, is.numeric, log)
    out2 <- gather(x2, variant2, `Cumulative Incidence`, -date)
    out2 <- select(out2, -date)

    out <- cbind(out, out2)
    gather(
      out, var, val, -date, -variant, -variant2
    )
  }
)


incid_plots <- map2(
  tall_incid, xaxis_breaks, function(x, breaks) {
    variants <- intersect(names(palette), unique(x$variant))
    values <- palette[variants]
    coeff <- max(x$val)

    ggplot(x) +
      geom_line(
        aes(date, val, col = variant, linetype = var),
        size = 1.1
      ) +
      scale_color_manual(
        values = values, labels = variant_nicenames
      ) +
      scale_x_date(
        breaks = breaks,
        date_labels = date_labels
      ) +
      scale_y_log10(
        labels = scales::label_log(digits = 1),
        sec.axis = dup_axis(
          name = "Cumulative Incidence",
          ##labels = trans_format('log10', math_format(10^.x))
        )
      ) +
      expand_limits(x = range(breaks)) +
      coord_cartesian(clip = "off") +
      ylab("Daily/Cumulative Incidence") +
      xlab("") +
  theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        ##Axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.box = "vertical",
        legend.margin = margin(),
        ##legend.position = c(0.2, 0.85),
        ## We don't actually want to show the right y-axis
        axis.line.y.right = element_line(color = "white"),
        axis.title.y.right = element_text(color = "white"),
        axis.text.y.right = element_text(color = "white"),
        axis.ticks.y.right = element_line(color = "white")
      )
  }
)

iwalk(
  incid_plots, function(p, name) {
    save_multiple(
      p, glue("figures/logscale_{name}_incidence")
    )
  }
)
