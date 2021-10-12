## orderly::orderly_develop_start()
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

### Panel B: Rt heatplots
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
pooled_estimate <- map2(
  epiestim_rt, select_variant, function(x, variants) {

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

rf <- colorRampPalette(rev(brewer.pal(8,'Spectral')))
r <- rf(32)

twodbin <- map2(
  pooled_estimate, select_variant,
  function(x, name) {
    xname <- glue(
      "Reproduction number for {to_title_case(name[1])}"
    )
    yname <- glue(
      "Reproduction number for {to_title_case(name[2])}"
    )
    maxrt <- ceiling(max(x, na.rm = TRUE))
    ggplot(x, aes(reference, variant)) +
      stat_bin2d(bins = 25) +
      scale_fill_gradientn(colours = r) +
      geom_abline(
        intercept = 0, slope = 1, col = 'grey50',
        linetype = "dashed", size = 1.2
      ) +
      xlim(0, maxrt) + ylim(0, maxrt) +
      xlab(xname) + ylab(yname) +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.key.width=unit(1.5,"cm")
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
