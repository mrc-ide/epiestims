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
