## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")
whole_cnty_time <- readRDS("epsilon_qntls_whole_country.rds")
x <- whole_cnty_time[["french"]]
whole_cnty_time[["french_betagamma"]] <- x[x$variant == "beta-gamma_vs_wild", ]
whole_cnty_time[["french"]] <- x[x$variant != "beta-gamma_vs_wild", ]
## Read all data first to get consistent ylimits
incidence <- readRDS("cuml_incid_all_variants.rds")
## Repeat Frnech data once to get betagamma as well.
incidence[["french_betagamma"]] <- incidence[["french"]]


## Epsilon estimates with non-overlapping weekly windows.
## Axis 2. Proportion of variant in the window used for estimation
nonovl_est <- readRDS("nonovl_epsilon_qntls.rds")
nonovl_est[["french_betagamma"]] <-
  nonovl_est[["french"]][nonovl_est[["french"]]$variant != "alpha_vs_wild", ]
nonovl_est[["french"]] <-
  nonovl_est[["french"]][nonovl_est[["french"]]$variant == "alpha_vs_wild", ]


## Epsilon estimates with non-overlapping weekly windows.
## and underreporting
nonovl_est_underrep <- readRDS("nonovl_epsilon_qntls_underrep.rds")
nonovl_est_underrep[["french_betagamma"]] <-
  nonovl_est_underrep[["french"]][nonovl_est[["french"]]$variant != "alpha_vs_wild", ]
nonovl_est_underrep[["french"]] <-
  nonovl_est_underrep[["french"]][nonovl_est_underrep[["french"]]$variant == "alpha_vs_wild", ]

## Put them together
nonovl_est <- map2(
  nonovl_est, nonovl_est_underrep, function(x, y) {
    bind_rows(None = x, `50%` = y, .id = "Underreporting")
  }
)



nonovl_prop <- readRDS("nonovl_prop_variant.rds")
## Has all the columns we need.
cols <- c("date", "wildtype", "betagamma",
          "cumulative_betagamma", "cumulative_wildtype",
          "proportion_betagamma")

nonovl_prop[["french_betagamma"]] <- map(
  nonovl_prop[["french"]], ~ .[, cols]
)

cols <- c("date", "wildtype", "alpha",
          "cumulative_alpha", "cumulative_wildtype",
          "proportion_alpha")

nonovl_prop[["french"]] <- map(nonovl_prop[["french"]], ~ .[, cols])

## For each variant, for each tmax we have the cumulative propotion
## in that week. Extract the last row.
nonovl_prop <- map(
  nonovl_prop, function(x) map_dfr(x, ~ .[nrow(.), ])
)

nonovl_prop <- map2(
  nonovl_prop,
  list(
    french = c("wildtype", "alpha"),
    uk_alpha_wild = c("cumulative_wildtype", "cumulative_alpha"),
    uk_delta_alpha = c("cumulative_alpha", "cumulative_delta"),
    french_betagamma = c("cumulative_wildtype", "cumulative_betagamma")
  ),
  function(x, cols) {
    out <- binconf(x[[cols[2]]], (x[[cols[1]]] + x[[cols[2]]]))
    cbind(x, data.frame(out))
  }
)


##
## broom::tidy(x)
plots2axis <- pmap(
  list(
    x = nonovl_est,
    y = nonovl_prop,
    w = whole_cnty_time,
    ## column = c("proportion_alpha", "proportion_alpha",
    ##             "proportion_delta", "proportion_betagamma"),
    y2label = c("Cumulative proportion of Alpha",
                "Cumulative proportion of Alpha",
                "Cumulative proportion of Delta",
                "Cumulative proportion of Beta/Gamma"),
    xmin = xaxis_breaks,
    col = palette[c("alpha", "alpha", "delta", "betagamma")]
  ),
  function(x, y, w, y2label, xmin, col) {
    message(y2label)
    x$date <- as.Date(x$date)
    z <- left_join(x, y, by = "date")
    w$date <- max(x$date) + 5 ## slightly to the right to avoid overlap
    ## coeff <-  max(z$`97.5%`) / max(z$proportion)
    coeff <- 1.8 ## For everyhing except delta
    if (x$variant[1] == "delta_vs_alpha") coeff <- 2.5
    z$proportion_scaled <- z$`PointEst` * coeff
    z$low_scaled <- z$`Lower` * coeff
    z$high_scaled <- z$`Upper` * coeff
    message("Max z$`97.5%` ", max(z$`97.5%`))
    message("Range of proportion ", range(z$proportion))
    message("Coeff = ", coeff)
    ggplot(z, aes(x = date)) +
      geom_point(
        aes(y = `50%`, shape = Underreporting), colour = col, size = 2,
        position = position_dodge(width = 2.5)
      ) +
      geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`, linetype = Underreporting),
        position = position_dodge(width = 2.5),
        size = 1.1, colour = col
      ) +
      geom_hline(
        yintercept = 1, linetype = "dashed", color = "red", size = 1.1
      ) +
      geom_point(aes(y = proportion_scaled), color = "blue") +
      geom_linerange(aes(ymin = low_scaled, ymax = high_scaled), color = "blue") +
      geom_line(aes(y = proportion_scaled), color = "blue", alpha = 0.5) +
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
      ## Add estimate for the entire country over the whole time period
      geom_point(
        data = w, aes(x = date, y = `50%`), fill = col, size = 3,
        shape = 23, col = "black", stroke = 2
      ) +
      geom_linerange(
        data = w, aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
        size = 1.1, colour = col
      ) +
      #########
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
    save_multiple(p, glue("figures/underrep_nonovl_{name}_2axis"))
  }
)













## Alpha multiplicative advantage by Volz et al
volzetal <- data.frame(
  date = as.Date("2020-12-31"),
  ymin = 1.4, ymax = 1.8, y = 1.74
)

##eps_over_time <- readRDS("epsilon_qntls_over_time.rds")
eps_over_time <- readRDS("epsilon_qntls_over_time.rds")
eps_over_time[["french_betagamma"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant != "alpha_vs_wild", ]
eps_over_time[["french"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant == "alpha_vs_wild", ]


rolling_window_prop <- readRDS("cuml_incid_all_variants.rds")
## Has all the columns we need.
cols <- c("date", "wildtype", "betagamma",
          "cumulative_betagamma", "cumulative_wildtype",
          "proportion_betagamma")

rolling_window_prop[["french_betagamma"]] <-
  rolling_window_prop[["french"]][, cols]


cols <- c("date", "wildtype", "alpha",
          "cumulative_alpha", "cumulative_wildtype",
          "proportion_alpha")

rolling_window_prop[["french"]] <-
  rolling_window_prop[["french"]][, cols]



rolling_window_prop <- map2(
  rolling_window_prop,
  list(
    french = c("wildtype", "alpha"),
    uk_alpha_wild = c("cumulative_wildtype", "cumulative_alpha"),
    uk_delta_alpha = c("cumulative_alpha", "cumulative_delta"),
    french_betagamma = c("cumulative_wildtype", "cumulative_betagamma")
  ),
  function(x, cols) {
    out <- binconf(x[[cols[2]]], (x[[cols[1]]] + x[[cols[2]]]))
    cbind(x, data.frame(out))
  }
)



######################################################################
## Plot with 2 y-axis
plots2axis <- pmap(
  list(
    x = eps_over_time,
    y = rolling_window_prop,
    ## column = c("proportion_alpha", "proportion_alpha",
    ##             "proportion_delta", "proportion_betagamma"),
    y2label =  c("Cumulative proportion of Alpha",
                 "Cumulative proportion of Alpha",
                 "Cumulative proportion of Delta",
                 "Cumulative proportion of Beta/Gamma"),
    xmin = xaxis_breaks,
    col = palette[c("alpha", "alpha", "delta", "betagamma")]
  ),
  function(x, y, y2label, xmin, col) {
    message(y2label)
    x$date <- as.Date(x$date)
    ##y <- select(y, date, PointEst:Upper)
    z <- left_join(x, y, by = "date")
    ## coeff <-  max(z$`97.5%`) / max(z$proportion)
    coeff <- 1.8 ## For everyhing except delta
    if (x$variant[1] == "delta_vs_alpha") coeff <- 2.5
    z$proportion_scaled <- z$`PointEst` * coeff
    z$low_scaled <- z$`Lower` * coeff
    z$high_scaled <- z$`Upper` * coeff
    message("Max z$`97.5%` ", max(z$`97.5%`))
    message("Range of proportion ", range(z$proportion))
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
      geom_point(aes(y = proportion_scaled), color = "blue") +
      geom_linerange(aes(ymin = low_scaled, ymax = high_scaled), color = "blue") +
      geom_line(aes(y = proportion_scaled), color = "blue", alpha = 0.5) +
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
    save_multiple(p, glue("figures/underrep_{name}_2axis"))
  }
)



