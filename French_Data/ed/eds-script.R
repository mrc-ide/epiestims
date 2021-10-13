library(tidyr)
library(tidyverse)
library(EpiEstim)

## read in england data
england <- readRDS("s_by_region_over25_pillar2_pcr.rds")
england <- england %>%
  select(specimen_date, nhser_name,
         s_positive_adj1, s_negative_adj1, s_na_adj1) %>%
  group_by(specimen_date, nhser_name) %>%
  summarise(s_positive_adj1 = sum(s_positive_adj1),
            s_negative_adj1 = sum(s_negative_adj1),
            s_na_adj1 = sum(s_na_adj1)) %>%
  mutate(region = tolower(gsub(" ", "_", nhser_name))) %>%
  rename(date = specimen_date)

## there are some samples not identified as negative or positive, we split those
## according to the ratio of the two. In the hopefully rare case that none are
## identified as negative or positive, these are split 50/50 - probably another
## way could be better?
na_positive_prop <- england$s_positive_adj1 /
  (england$s_positive_adj1 + england$s_negative_adj1)
na_positive_prop <- ifelse(is.na(na_positive_prop), 0.5, na_positive_prop)
s_positive_adj1_extra <- round(na_positive_prop * england$s_na_adj1)
s_negative_adj1_extra <- england$s_na_adj1 - s_positive_adj1_extra
england$s_positive_adj1 <- england$s_positive_adj1 + s_positive_adj1_extra
england$s_negative_adj1 <- england$s_negative_adj1 + s_negative_adj1_extra


## Run England fits in two phases:
## 1. Positive corresponds to wildtype, Negative to B.1.1.7
## 2. Negative corresponds to B.1.1.7, Positive to B.1.617.2

inputs_phase_1 <- list(fit1 = list(date_start = "2020-09-01",
                                   date_end = "2021-02-01",
                                   t_min = NULL),
                       fit2 = list(date_start = "2020-11-01",
                                   date_end = "2021-01-01",
                                   t_min = NULL),
                       fit3 = list(date_start = "2020-11-01",
                                   date_end = "2021-01-01",
                                   t_min = 14L))

england_phase_1 <-
  lapply(inputs_phase_1,
         function(input) {
           epiestim_fit(england, sircovid::regions("england"),
                        input$date_start, input$date_end, input$t_min,
                        "s_positive_adj1", "s_negative_adj1")})

inputs_phase_2 <- list(fit1 = list(date_start = "2021-04-01",
                                   date_end = max(england$date),
                                   t_min = NULL),
                       fit2 = list(date_start = "2021-05-01",
                                   date_end = max(england$date),
                                   t_min = NULL),
                       fit3 = list(date_start = "2021-04-01",
                                   date_end = max(england$date),
                                   t_min = 14L),
                       fit4 = list(date_start = "2021-04-15",
                                   date_end = max(england$date),
                                   t_min = 14L))

england_phase_2 <-
  lapply(inputs_phase_2,
         function(input) {
           epiestim_fit(england, sircovid::regions("england"),
                        input$date_start, input$date_end, input$t_min,
                        "s_negative_adj1", "s_positive_adj1")})


dir.create("outputs")

saveRDS(england_phase_1, "outputs/england_phase_1.rds")
saveRDS(england_phase_2, "outputs/england_phase_2.rds")

png(filename = "outputs/england_prop_positive_phase_1.png",
    width = 1200, height = 600)
plot_england_prop_positive(england, "2020-09-01", "2021-03-15")
dev.off()

png(filename = "outputs/england_prop_positive_phase_2.png",
    width = 1200, height = 600)
plot_england_prop_positive(england, "2021-03-15", max(england$date))
dev.off()


# ## read in french data
# france <- read.csv("rtm_france_wide.csv")
# france <- france %>%
#   rename(region = reg_name) %>%
#   mutate(n_wildtype = round(france$n_tests_no_variant_dep_7day_allage / 7)) %>%
#   mutate(n_uk = round(france$n_tests_uk_susp_dep_7day_allage / 7)) %>%
#   mutate(n_sa_brazil = round(france$n_tests_sa_brazil_susp_dep_7day / 7)) %>%
#   filter(!is.na(n_wildtype))
#
#
# ## French regions
# metropolitan_france <- c("Ile-de-France", "Centre-Val de Loire",
#                          "Bourgogne-Franche-Comte", "Normandy",
#                          "Hauts-de-France", "Grand Est", "Pays de la Loire",
#                          "Brittany", "Nouvelle-Aquitaine", "Occitanie",
#                          "Auvergne-Rhone-Alpes", "Provence-Alpes-Cote d'Azur",
#                          "Corsica")
# all_france <- c(metropolitan_france, "Guadeloupe", "Martinique",
#                 "French Guiana", "Reunion", "Mayotte")
#
# ## Just fit to metropolitan France
# france_metro_fit <- epiestim_fit(france, metropolitan_france,
#                                  "2021-02-18", max(france$date),
#                                  "n_wildtype", "n_uk", "n_sa_brazil")
#
# ## fit to overseas regions as well
# france_all_fit <- epiestim_fit(france, all_france,
#                                "2021-02-18", max(france$date),
#                                "n_wildtype", "n_uk", "n_sa_brazil")
#
# png(filename = "outputs/france_metro_fit.png", width = 960, height = 960)
# plot_france(france_metro_fit)
# dev.off()
#
# png(filename = "outputs/france_all_fit.png", width = 960, height = 960)
# plot_france(france_all_fit)
# dev.off()
#
# saveRDS(france_metro_fit, "outputs/france_metro_fit.rds")
# saveRDS(france_all_fit, "outputs/france_all_fit.rds")
