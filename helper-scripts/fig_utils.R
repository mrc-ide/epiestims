scenario_type_labeller <- function(x) {
  x$scenario_type <- case_when(
    x$label == "X 0.5" ~ "Low",
    x$label == "X 1.5" ~ "Moderate",
    x$label == "X 1" ~ "Baseline",
    x$label == "X 2" ~ "High",
    ## Vary offspring
    x$label == "0.1" ~ "High",
    x$label == "0.5" ~ "Moderate",
    x$label == "1" ~ "Low",
    ## Underreporting
    x$label == "p0.2" ~ "0.2",
    x$label == "p0.5" ~ "0.5",
    x$label == "p0.8" ~ "0.8"
  )
  x
}



multiplier_label <- function(val, ref) {
  paste("X", round(val/ref, 1))
}

theme_manuscript <- function(base_size = 16) {
  theme_classic() %+replace%
    theme(
      text = element_text(size = base_size),
      ##axis.title = element_text(size = 12),
      ##axis.line = element_line(size = 1.05),
      ##strip.background = element_rect(size = 1.05),
      legend.position = "top",
      axis.text.x = element_text(angle = 55)
    )
}
## give filename without the extension
save_multiple <- function(plot, filename) {
  ggsave(
    filename = glue("{filename}.pdf"),
    plot
  )
  ggsave(
    filename = glue("{filename}.png"),
    plot)
}

rt_labeller <- function(val) {
  paste("Reference Rt:", val)
}

rt_change_labeller <- function(val) {
  paste("Reference Rt change:", val)
}

tmax_labeller <- function(val) {
  paste(val, "days")
}

## df is the output of summarise_median_err
format_median_err <- function(df) {
  df$formatted <-
  glue("{df$median_med}",
       " ({df$median_low}, {df$median_high})")
  df <- select(df, true_eps, tmax, formatted) %>%
    spread(key = tmax, value = formatted)
  df
}

pretty_ci <- function(val, low, high, round_to = 2) {
  f <- function(x) {
    format(round(x, round_to), nsmall = 2)
  }
  glue("{f(val)} \n ({f(low)}, {f(high)})")
}

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
