## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")
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
### Panel B: Rt heatplots
#################################################
infiles <- list(
  "Rt_epi_fr.rds", "Rt_epi_UK1.rds", "Rt_epi_UK2.rds",
  "Rt_epi_fr.rds"
)
names(infiles) <- c(
  "french", "uk_alpha_wild", "uk_delta_alpha", "french_betagamma"
)
epiestim_rt <- map(infiles, readRDS)
## Each element of epiestim_rt is a list with
## elements corresponding to variants.
## Here we prepare a list to select the appropriate
## variant outputs
select_variant <- list(
  c("wild", "alpha"), c("wild", "alpha"),
  c("alpha", "delta"), c("wild", "beta/gamma")
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
    if (name[2] == "beta/gamma") yname <- "Reproduction number for Beta/Gamma"
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
    xlim <- 3
    if (name[2] == "delta") xlim <- 6
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
      xlim(0, xlim) +
      ylim(0, xlim) +
      xlab(xname) +
      ylab(yname) +
     coord_fixed() +
      labs("Density") +
      theme_manuscript() +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.key.height = unit(1.2, "cm"),
        legend.position = "right"
      )
  }
)

iwalk(
  twodbin, function(p, name) {
    save_multiple(
      p, glue("figures/{name}_2dbin")
    )
    ##knitr::plot_crop(glue("figures/{name}_2dbin.png"))
  }
)



