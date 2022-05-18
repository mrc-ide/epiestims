## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")
regional <- readRDS("epsilon_qntls_per_region.rds")
regional[["french_betagamma"]] <-
  regional[["french"]][regional[["french"]]$variant != "alpha_vs_wild", ]
regional[["french"]] <-
  regional[["french"]][regional[["french"]]$variant == "alpha_vs_wild", ]

national <- readRDS("epsilon_qntls_whole_country.rds")
national[["french_betagamma"]] <-
  national[["french"]][national[["french"]]$variant != "alpha_vs_wild", ]
national[["french"]] <-
  national[["french"]][national[["french"]]$variant == "alpha_vs_wild", ]

national <- map2(
  national, list("France", "England", "England", "France"),
  function(x, y) {
    x$region <- y
    x
  }
)


#################################################
###### Panel C. Regional estimates
#################################################
regional_plots <- pmap(
  list(x = regional, y = national,
       z = palette[c("alpha", "alpha", "delta", "betagamma")]
       ),
  function(x, y, z) {
    ## Easier than creating a palette
    x$shape <- 21
    y$shape <- 23
    x$colour <- z
    y$colour <- "black"
    ## Arrange columns in the same order
    ## for rbind
    y <- y[, colnames(x)]
    ## x <- arrange(x, desc(`50%`))
    x <- arrange(x, `50%`)
    x <- rbind(y, x)
    x$region <- factor(
      x$region,
      levels = x$region,
      ordered = TRUE
    )
    xmax <- 1.8 ## For everyhing except delta
    if (x$variant[1] == "delta_vs_alpha") xmax <- 2.5
    ggplot(x) +
      geom_linerange(
        aes(region, ymin = `2.5%`, ymax = `97.5%`),
        size = 1.1, colour = z
      ) +
      geom_point(
        aes(region, `50%`, shape = shape, colour = colour),
        size = 2, fill = z, stroke = 2
      ) +
      geom_hline(
        yintercept = 1, linetype = "dashed", color = "red",
        size = 1.2
      ) +
      ##expand_limits(y = 1) +
      ##ylim(0.5, 2) +
      scale_shape_identity(
        breaks = c(21, 23)
      ) +
      scale_colour_identity(
        breaks = c(z, "black")
      ) +
      scale_x_discrete(
        labels = region_short_names
      ) +
      scale_y_continuous(
        limits = c(0, xmax),
        breaks = seq(0, xmax, by = 0.5)
      ) +
      ylab("Effective transmission advantage") +
      ##coord_flip() +
      theme_manuscript() +
      theme(
        ##axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()
      )
  }
)

iwalk(
  regional_plots, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_regional")
    )
  }
)

#################################
## Non-overlapping regional estimates
#################################
infiles <- list(
  french = 'I_fr.rds', uk_alpha_wild = 'I_UK1.rds',
  uk_delta_alpha = 'I_UK2.rds'
)

incidence <- map(infiles, readRDS)
nonovl <- readRDS("nonovl_weekly_regional_epsilon_qntls.rds")

nonovl <- map2(
  incidence, nonovl, function(x, y) {
    incid <- x[[1]]
    y$date <- incid$date[as.integer(y$tmax)]
    y
  }
)
x <- nonovl[["french"]]
nonovl[["french_betagamma"]] <- x[x$variant == "beta-gamma_vs_wild", ]
nonovl[["french"]] <- x[x$variant != "beta-gamma_vs_wild", ]
palette2 <- palette[c("alpha", "alpha", "delta", "betagamma")]

plots <- pmap(
  list(
    x = nonovl, row = c(3, 2, 2, 3), z = palette2, breaks = xaxis_breaks
  ),
  function(x, row, z, breaks) {
  x$date <- as.Date(x$date)
  x$location_short <- region_short_names(x$location)
  p <- ggplot(x) +
    geom_linerange(
      aes(
        x = date, ymin = `2.5%`, ymax = `97.5%`
      ),
        size = 1.1, colour = z
    ) +
    geom_point(
      aes(x = date, y = `50%`),
      colour = z
    )  +
    geom_hline(
      yintercept = 1, linetype = "dashed", color = "red",
      size = 1.2
    ) +
    facet_wrap(
      ~location_short, nrow = row, strip.position = "top"
    ) +
    scale_x_date(breaks = breaks, date_labels = date_labels) +
    ylab("Effective transmission advantage") +
    theme_manuscript() +
    theme(
      axis.title.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      ##ggh4x.axis.nesttext.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  p

})



iwalk(
  plots, function(p, name) {
    save_multiple(
      p, glue("figures/nonovl_{name}_regional")
    )
  }
)
