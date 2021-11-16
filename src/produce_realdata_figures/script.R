## orderly::orderly_develop_start()
source("R/fig_utils.R")
dir.create("figures")
mypercent <- function(x) scales::percent(x)
palette <- c(
  wildtype = "#0f0e0e",
  alpha = "#56B4E9",
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
  delta = "Delta", betagamma = "Beta/Gamma",
  `beta/gamma` = "Beta/Gamma"
)

## Read all data first to get consistent ylimits
incidence <- readRDS("cuml_incid_all_variants.rds")
## Repeat Frnech data once to get betagamma as well.
incidence[["french_betagamma"]] <- incidence[["french"]]
date_min <- map(incidence, ~ min(.$date))
#### Date breaks to start in the 1st of every month
xaxis_breaks <- map(
  incidence, function(x) {
    xmin <- round_date(min(x$date), "14 days")
    xmax <- round_date(max(x$date), "month")
    seq(as.Date(xmin), as.Date(xmax), "2 months")
  }
)

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
  short_name <- case_when(
    region == "Auvergne-Rhône-Alpes" ~ "ARA",
    region ==  "Bourgogne-Franche-Comté" ~ "BFC",
    region ==  "Bretagne"  ~ "BRE",
    region ==  "Centre-Val de Loire" ~ "CVL",
    region ==  "Corse" ~ "20R",
    region ==  "Grand Est" ~ "GES",
    region ==  "Guadeloupe" ~ "GP",
    region ==  "Guyane"  ~ "GF",
    region ==  "Hauts-de-France" ~ "HDF",
    region ==    "Île-de-France" ~ "IDF",
    region ==     "La Réunion" ~ "RE",
    region ==    "Martinique" ~ "MQ",
    region ==    "Mayotte" ~ "YT",
    region ==    "Normandie" ~ "NOR",
    region ==    "Nouvelle-Aquitaine" ~ "NAQ",
    region ==    "Occitanie" ~ "OCC",
    region ==    "Pays de la Loire" ~ "PDL",
    region ==    "Provence-Alpes-Côte d'Azur" ~ "PAC",
    region ==    "East of England" ~ "EE",
    region ==    "London" ~ "LON",
    region ==    "Midlands" ~ "MID",
    region ==    "North East and Yorkshire" ~ "NE",
    region ==    "North West" ~ "NW",
    region ==    "South East" ~ "SE",
    region ==    "South West" ~ "SW",
    region ==    "England" ~ "ENG",
    region ==    "France" ~ "FR"
  )
  short_name
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
regional_incid[["england"]] <- rbind(
  regional_incid[["uk1"]], regional_incid[["uk2"]]
)

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
    ## x <- arrange(x, desc(`50%`))
    x <- arrange(x, `50%`)
    x <- rbind(y, x)
    x$region <- factor(
      x$region,
      levels = x$region,
      ordered = TRUE
    )
    xmax <- 1.8 ## For everyhing except delta
    if (x$variant[1] == "delta_vs_alpha") xmax <- 2
    ggplot(x) +
      geom_point(
        aes(region, `50%`, color = color),
        size = 4
      ) +
      geom_linerange(
        aes(region, ymin = `2.5%`, ymax = `97.5%`, color = color),
        size = 1.1
      ) +
      geom_hline(
        yintercept = 1, linetype = "dashed", color = "red",
        size = 1.2
      ) +
      ##expand_limits(y = 1) +
      ##ylim(0.5, 2) +
      scale_color_identity(
        breaks = c("#0f0e0e", palette[[y$region[1]]])
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
    ggplot(x) +
      geom_bin2d(
        aes(reference, variant, fill = ..density..),
        bins = 25) +
      scale_fill_gradientn(colours = r) +
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
      ##scale_fill_gradientn(colours = r) +
      xlim(0, 3) +
      ylim(0, 3) +
      xlab(xname) +
      ylab(yname) +
      labs("Density") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.key.width = unit(2.5, "cm")
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
  eps_over_time, xaxis_breaks, function(x, xmin) {
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
        breaks = xmin,
        ## date_breaks = date_breaks,
        date_labels = date_labels
        ## limits = c(as.Date(xmin), NA)
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
        ylab("Effective transmission advantage") +
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
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1.2) +
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
    if (x$variant[1] == "delta_vs_alpha") coeff <- 2

    z$proportion_scaled <- z$proportion * coeff
    message("Max z$`97.5%`", max(z$`97.5%`))
    message("Range of proportion", range(z$proportion))
    message("Coeff = ", coeff)
    ggplot(z, aes(x = date)) +
      geom_point(
        aes(y = `50%`), size = 2, colour = col
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
        limits = c(0, 1.8),
        breaks = seq(0, 1.8, by = 0.5)
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

### Summary tables
naive_eps <- map(
  c(french = "naive_epsilon_fr.rds",
    uk_alpha_wild = "naive_epsilon_UK1.rds",
    uk_delta_alpha = "naive_epsilon_UK2.rds"), function(x) {
      readRDS(x)$sum_epsilon
    }
)

naive_eps[["french_betagamma"]] <- naive_eps[["french"]][["wild_beta/gamma"]]
naive_eps[["french"]] <- naive_eps[["french"]][["wild_alpha"]]
naive_eps[["uk_alpha_wild"]] <- naive_eps[["uk_alpha_wild"]][[1]]
naive_eps[["uk_delta_alpha"]] <- naive_eps[["uk_delta_alpha"]][[1]]

## MV-EpiEstim time periods Q1-4
mvepi_q <- readRDS("epsilon_qntls_time_periods.rds")
mvepi_q[["french_betagamma"]] <-
  mvepi_q[["french"]][mvepi_q[["french"]]$variant != "alpha_vs_wild", ]
mvepi_q[["french"]] <-
  mvepi_q[["french"]][mvepi_q[["french"]]$variant == "alpha_vs_wild", ]

pretty_ci <- function(val, low, high, round_to = 2) {
  f <- function(x) {
    format(round(x, round_to), nsmall = 2)
  }
  glue("{f(val)} ({f(low)}, {f(high)})")
}

## get dates corresponding to quarters
periods <- readRDS("periods.rds")
periods[["french_betagamma"]] <- periods[["periods_fr"]]

quarter_dates <- map2(
  incidence, periods, function(incid, period) {
    t_min <- period$intervals[-1] ## Remove the first element
    t_min <- head(t_min, -1)
    t_max <- period$intervals[-c(1, 2)]
    t_min <- format(incid$date[t_min], "%d-%b-%Y")
    t_max <- format(incid$date[t_max], "%d-%b-%Y")
    glue("{t_min} to {t_max}")
  }
)


si_tables <- pmap(
  list(naive = naive_eps, mvepi = regional, overall = national, qx = mvepi_q,
       qd = quarter_dates),
  function(naive, mvepi, overall, qx, qd) {
    overall$region <- "All"
    mvepi <- rbind(overall[, colnames(mvepi)], mvepi)
    qx$region <- qx$time_period
    mvepi <- rbind(qx[, colnames(mvepi)], mvepi)
    idx <- which(is.na(naive$med))
    naive$formatted <- pretty_ci(naive$med, naive$low, naive$up)
    naive$formatted[idx] <- "-"
    mvepi$formatted <- pretty_ci(mvepi$`50%`, mvepi$`2.5%`, mvepi$`97.5%`)
    naive <- select(naive, "Region/Time-Period" = name, `Naive` = formatted)
    mvepi <- select(mvepi, "Region/Time-Period" = region, `MV-EpiEstim` = formatted)
    out <- left_join(naive, mvepi, by = "Region/Time-Period")
    out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 1"] <- qd[1]
    out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 2"] <- qd[2]
    out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 3"] <- qd[3]
    out["Region/Time-Period"][out["Region/Time-Period"] == "Quarter 4"] <- qd[4]
    out
  }
)


iwalk(
  si_tables, function(x, i) {
    out <- ggtexttable(
      x, rows = NULL, cols = colnames(x),
      theme = ttheme("mBlue", base_size = 9)
    )
    save_multiple(out, glue("figures/{i}_si_table"))
  }
)
