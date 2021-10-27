## orderly::orderly_develop_start()
source("R/fig_utils.R")
dir.create("figures")

palette <- c(
  wildtype = "#0f0e0e",
  alpha = "#E69F00",
  betagamma = "#56B4E9",
  `beta/gamma` = "#56B4E9",
  delta = "#009E73",
  England = "#cc0000",
  France = "#0000ff"
)

date_breaks <- "6 weeks"
date_labels <- "%d-%b"


variant_nicenames <- c(
  wildtype = "Wildtype", alpha = "Alpha",
  delta = "Delta", betagamma = "Beta/Gamma"
)

## Read all data first to get consistent ylimits
incidence <- readRDS("cuml_incid_all_variants.rds")
## Repeat Frnech data once to get betagamma as well.
incidence[["french_betagamma"]] <- incidence[["french"]]
date_min <- map(incidence, ~ min(.$date))

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

## Alpha multiplicative advantage by Volz et al
volzetal <- data.frame(
  date = as.Date("2020-12-31"),
  ymin = 1.4, ymax = 1.8, y = 1.74
)

eps_over_time <- readRDS("epsilon_qntls_over_time.rds")
eps_over_time[["french_betagamma"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant != "alpha_vs_wild", ]
eps_over_time[["french"]] <-
  eps_over_time[["french"]][eps_over_time[["french"]]$variant == "alpha_vs_wild", ]

## Region short names
region_short_names <- function(region) {
    lookup <- c(ARA = "Auvergne-Rhône-Alpes",
                BFC = "Bourgogne-Franche-Comté",
                BRE = "Bretagne",
      CVL = "Centre-Val de Loire",
      `20R` = "Corse",
      GES = "Grand Est",
      GP = "Guadeloupe",
      GF = "Guyane",
      HDF = "Hauts-de-France",
      IDF = "Île-de-France",
      RE  = "La Réunion",
      MQ = "Martinique",
      YT = "Mayotte",
      NOR = "Normandie",
      NAQ = "Nouvelle-Aquitaine",
      OCC = "Occitanie",
      PDL = "Pays de la Loire",
      PAC = "Provence-Alpes-Côte d'Azur",
      EE = "East of England",
      LON = "London",
      MID = "Midlands",
      NE = "North East and Yorkshire",
      NW = "North West",
      SE = "South East",
      SW = "South West",
      ENG = "England",
      FR = "France"
      )
    names(lookup)[lookup %in% region]
}


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
        date_breaks = date_breaks,
        date_labels = date_labels
      ) +
      coord_cartesian(clip = "off") +
      ylab("Daily incidence") +
      xlab("") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        ##axis.title.x = element_blank(),
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
#################################################
###### Panel C. Regional estimates
#################################################
regional_plots <- map2(
  regional, national, function(x, y) {
    ## Easier than creating a palette
    x$color <- "#0f0e0e" ## muted black
    y$color <- palette[[y$region[1]]]
    ## Arrange columns in the same order
    ## for rbind
    y <- y[, colnames(x)]
    x <- arrange(x, desc(`50%`))
    x <- rbind(x, y)
    x$region <- factor(
      x$region,
      levels = x$region,
      ordered = TRUE
    )

    ggplot(x) +
      geom_point(
        aes(region, `50%`, color = color),
        size = 4
      ) +
      geom_linerange(
        aes(region, ymin = `2.5%`, ymax = `97.5%`, color = color),
        size = 1.1
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      ##expand_limits(y = 1) +
      ##ylim(0.5, 2) +
      scale_color_identity(
        breaks = c("#0f0e0e", palette[[y$region[1]]])
      ) +
      scale_x_discrete(
        labels = region_short_names
      ) +
      scale_y_continuous(
        limits = c(0, 2),
        breaks = seq(0, 2, by = 0.5)
      ) +
      ylab("Effective transmission advantage") +
      ##coord_flip() +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
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
#################################################
### Panel B: Rt heatplots
#################################################
infiles <- list(
  "Rt_epi_fr.rds", "Rt_epi_UK1.rds", "Rt_epi_UK2.rds",
  "Rt_epi_fr.rds"
)
names(infiles) <- names(incidence)
epiestim_rt <- map(infiles, readRDS)
## Each element of epiestim_rt is a list with
## elements corresponding to variants.
## Here we prepare a list to select the appropriate
## variant outputs
select_variant <- list(
  c("wild", "alpha"), c("wild", "alpha"),
  c("alpha", "delta"), c("alpha", "beta/gamma")
)
names(select_variant) <- names(epiestim_rt)
## pooled estimates from region-weeks included
pooled_estimate <- pmap(
  list(
    x = epiestim_rt, variants = select_variant
  ),
  function(x, variants) {
    reference <- x[[variants[[1]]]]
    variant <- x[[variants[[2]]]]
    ## Rt_incl is an indicator variable
    ## for when Rt estimate from a given day and
    ## region is  included. exclude if NA
    ## include if 1.
    remove <- reference$Rt_incl[, -1] * variant$Rt_incl[, -1]
    ## Rt samples. Dimensions T X Regions X samples
    ref_samples <- reference$Rt_s
    ## Remove the ones you don't want
    for (i in 1:dim(ref_samples)[3]) {
      ref_samples[, , i] <- as.matrix(ref_samples[, , i] * remove)
    }
    ## Similarly for variant
    var_samples <- variant$Rt_s
    ## Remove the ones you don't want
    for (i in 1:dim(var_samples)[3]) {
      var_samples[, , i] <- as.matrix(var_samples[, , i] * remove)
    }
    data.frame(
      reference = c(ref_samples), variant = c(var_samples)
    )
  }
)

rf <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
r <- rf(32)

## Plotting eptimates from MV-EpiEstim on the
## 2D bin plot
## x2 <- seq(0,as.numeric(max(x,na.rm=TRUE)),length.out=10)
## d2 <- data.frame(x = x2,
##                  y = as.numeric(epsilon[1])*x2,
##                  ylower = as.numeric(epsilon[2])*x2,
##                  yupper = as.numeric(epsilon[3])*x2)


twodbin <-pmap(
  list(
    x = pooled_estimate, name = select_variant,
    mv_estimate = national
  ), function(x, name, mv_estimate) {
    if (name[1] == "wild") name[1] <- "wildtype"
    xname <- glue(
      "Reproduction number for {to_title_case(name[1])}"
    )
    yname <- glue(
      "Reproduction number for {to_title_case(name[2])}"
    )
    message(yname)
    maxrt <- ceiling(max(x, na.rm = TRUE))
    x2 <- seq(0, maxrt, length.out = 10)
    mv_estim <- data.frame(
      x = x2,
      y = x2 * mv_estimate[["50%"]],
      ymin = x2 * mv_estimate[["2.5%"]],
      ymax = x2 * mv_estimate[["97.5%"]]
    )
    maxrt <- max(c(maxrt, mv_estim$ymax))
    ggplot() +
      geom_bin2d(
        data = x, aes(reference, variant), bins = 25
      ) +
      geom_line(
        data = mv_estim,
        aes(x = x, y = y),
        col = palette[[name[2]]]
      ) +
      geom_ribbon(
        data = mv_estim,
        aes(x = x, ymin = ymin, ymax = ymax),
        fill = palette[[name[2]]], alpha = 0.3
      ) +
      geom_abline(
        intercept = 0, slope = 1, col = "grey50",
        linetype = "dashed", size = 1.2
      ) +
      scale_fill_gradientn(colours = r) +
      xlim(0, maxrt) +
      ylim(0, maxrt) +
      xlab(xname) +
      ylab(yname) +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.key.width = unit(1.5, "cm")
      )
  }
)

iwalk(
  twodbin, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_2dbin")
    )
  }
)



################################################
###### Panel D. Estimates over time
################################################

plots_over_time <- map2(
  eps_over_time, date_min, function(x, xmin) {
    x$date <- as.Date(x$date)
    ggplot(x) +
      geom_point(aes(date, `50%`), size = 2) +
      geom_linerange(
        aes(date, ymin = `2.5%`, ymax = `97.5%`),
        size = 1.1
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_y_continuous(
        limits = c(0.5, 2),
        breaks = seq(0.5, 2, by = 0.5)
      ) +
      scale_x_date(
        date_breaks = date_breaks,
        date_labels = date_labels,
        limits = c(as.Date(xmin), NA)
      ) +
      ylab("Effective transmission advantage") +
      xlab("Estimation using data reported up to") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        ##axis.title.x = element_blank(),
        legend.title = element_blank()
      )
  }
)

iwalk(
  plots_over_time, function(p, name) {
    if (name == "uk_alpha_wild") {
      p <- p +
        geom_point(
          data = volzetal, aes(date, y),
          col = "#E69F00", size = 2,
          position = position_nudge(x = 2)
        ) +
        geom_linerange(
          data = volzetal,
          aes(date, ymin = ymin, ymax = ymax),
          col = "#E69F00", size = 1.1,
          position = position_nudge(x = 2)
      )
    }
    save_multiple(
      p, glue("figures/{name}_over_time")
    )
  }
)

##################################################
###### Panel E. Estimates with proportion of variant
##################################################
mypercent <- function(x) scales::percent(x, accuracy = 0.1)
eps_over_time_with_prop <- readRDS("epsilon_estimates_with_variant_proportion.rds")
eps_over_time_with_prop[["french_betagamma"]] <-
  eps_over_time_with_prop[["french"]][eps_over_time_with_prop[["french"]]$variant != "alpha_vs_wild", ]
eps_over_time_with_prop[["french"]] <-
  eps_over_time_with_prop[["french"]][eps_over_time_with_prop[["french"]]$variant == "alpha_vs_wild", ]

## custom_eps_over_time <- readRDS(
##   "custom_epsilon_estimates_with_variant_proportion.rds"
## )

## custom_priors <- readRDS("custom_epsilon_priors.rds")

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

eps_with_prop <- map2(
    both_together,
    c("Proportion of Alpha",
      "Proportion of Alpha",
      "Proportion of Delta",
      "Proportion of Beta/Gamma"
      ), function(x, xlabel) {
        message(xlabel)
        ##x <- x[x$prior == "Default priors", ]
        p <- ggplot(x) +
          geom_point(
            aes(proportion, `50%`), size = 2
          ) +
          geom_linerange(
            aes(proportion, ymin = `2.5%`, ymax = `97.5%`),
            size = 1.1
          ) +
          geom_hline(
            yintercept = 1, linetype = "dashed", color = "red"
          ) +
        expand_limits(y = 1) +
        scale_x_continuous(labels = mypercent) +
        ylab("Effective Transmission Advantage") +
        xlab(xlabel) +
        theme_manuscript() +
        theme(
          axis.text.x = element_text(angle = 0),
          legend.title = element_blank()
        )
    p
  }
)

iwalk(
  eps_with_prop, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_over_proportion")
    )
  }
)

## Somewhat tricky to see the different prior is
## making
eps_with_prop <- map2(
    both_together,
    c("Proportion of Alpha",
      "Proportion of Alpha",
      "Proportion of Delta",
      "Proportion of Beta/Gamma"
      ), function(x, xlabel) {
        message(xlabel)
        xmax <- min(c(max(x$proportion)/2, 0.1))
        ##x$proportion <- log(x$proportion)
        p <- ggplot(x) +
      geom_ribbon(
        aes(proportion, ymin = `2.5%`, ymax = `97.5%`),
        alpha = 0.2
      ) +
      geom_line(
        aes(proportion, `50%`), size = 1.1
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      expand_limits(y = 1) +
      scale_x_continuous(labels = mypercent) +
      ## scale_color_manual(
      ##   values = c(
      ##     `Default prior` = "#0f0e0e",
      ##     `Informative prior` = "#CC79A7"
      ##   ),
      ##   aesthetics = c("col", "fill")
      ## ) +
      ylab("Effective Transmission Advantage") +
      xlab(xlabel) +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.title = element_blank()
      )

    p
  }
)

iwalk(
  eps_with_prop, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_over_proportion_SI")
    )
  }
)
### Putting everthing together
## p1 <- incid_plots[[1]]
## p2 <- twodbin[[1]]
## p3 <- regional_plots[[1]]
## p4 <- plots_over_time[[1]]
## p5 <- eps_with_prop[[1]]

## p <- (p1 | p2 | p3) / (p4 | p5)
## save_multiple(p, "test.pdf")

######################################################################
## Plot with 2 y-axis
plots2axis <- pmap(
  list(
    x = eps_over_time,
    y = both_together,
    y2label =  c("Proportion of Alpha",
                 "Proportion of Alpha",
                 "Proportion of Delta",
                 "Proportion of Beta/Gamma"),
    xmin = date_min
  ),
  function(x, y, y2label, xmin) {
    y <- select(y, date, proportion)
    x$date <- as.Date(x$date)
    z <- left_join(x, y, by = "date")
    coeff <-  max(z$`97.5%`) / max(z$proportion)
    z$proportion_scaled <- z$proportion * coeff
    message("Max z$`97.5%`", max(z$`97.5%`))
    message("Range of proportion", range(z$proportion))
    message("Coeff = ", coeff)
    ggplot(z, aes(x = date)) +
      geom_point(
        aes(y = `50%`), size = 2
      ) +
      geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`),
        size = 1.1
      ) +
      geom_hline(
        yintercept = 1, linetype = "dashed", color = "red"
      ) +
      geom_line(aes(y = proportion_scaled), linetype = "dashed") +
      scale_y_continuous(
        sec.axis = sec_axis(~./coeff, name = y2label),
        limits = c(0, 2),
        breaks = seq(0, 2, by = 0.5)
      ) +
      scale_x_date(
        date_breaks = date_breaks,
        date_labels = date_labels,
        limits = c(as.Date(xmin), NA)
      ) +
      coord_cartesian(clip = "off") +
      ylab("Effective Transmission Advantage") +
      xlab("") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        legend.title = element_blank(),
        axis.line.y.right = element_line(linetype = "dashed")
      )

  }
)

iwalk(
  plots2axis, function(p, name) {
    if (name == "uk_alpha_wild") {
      p <- p +
        geom_point(
          data = volzetal, aes(date, y),
          col = "#E69F00", size = 2,
          position = position_nudge(x = 2)
        ) +
        geom_linerange(
          data = volzetal,
          aes(date, ymin = ymin, ymax = ymax),
          col = "#E69F00", size = 1.1,
          position = position_nudge(x = 2)
      )
    }
    save_multiple(p, glue("figures/{name}_2axis"))
  }
)

