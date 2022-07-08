library(orderly)
library(purrr)


###########################
## ANALYSIS OF REAL DATA ##
###########################

## Analysis, results and figures for main text (Fig 1) and
## SI Secs 2 & 4

a <- orderly_run("clean_french_england_data")
orderly_commit(a)

a <- orderly_run("naive_epsilon_estimates")
orderly_commit(a)

a <- orderly_run("mv_epiestim_weekly_estimates")
orderly_commit(a)

a <- orderly_run("src/produce_realdata_figures/")

a <- orderly_run("mv_epiestim_regional_estimates")
orderly_commit(a)

a <- orderly_run("src/produce_regional_eps_figures/")



################################
## ANALYSIS OF SIMULATED DATA ##
################################

## Analysis, results and figures for main text (Fig 2) and
## SI Sec 5

## Function running three orderly tasks for a given scenario
# Tasks are:
# 1) simulation of incidence data;
# 2) estimation of transmission advantage;
# 3) summary of MV-EpiEstim output

scenario_tasks <- function(scenario_name) {
  
  sim_task <- paste0("sims_", scenario_name)
  estimation_task <- scenario_name
  summary_task <- paste0("summarise_", scenario_name)
  
  a <- orderly_run(sim_task, parameters = list(short_run = FALSE))
  orderly_commit(a)
  
  a <- orderly_run(estimation_task, parameters = list(short_run = FALSE))
  orderly_commit(a)
  
  a <- orderly_run(summary_task)
  orderly_commit(a)
  
}


## Apply scenario_tasks function to scenarios featured in our analysis


# SI distribution of variant has different mean
# (See SI Sec 5.4)
scenario_tasks("one_location_vary_si")

# SI distribution of variant has different CV
# (See SI Sec 5.6)
scenario_tasks("one_location_vary_cv")

# Use a negative binomial offspring distribution
# (See SI Sec 5.8)
scenario_tasks("one_location_vary_offspring")

# Assess sensitivity of MV-EpiEstim to under-reporting
# (See SI Sec 5.9)
scenario_tasks("one_location_vary_underreporting")

# Time-varying Rt profiles in one or two locations
# (See SI Secs 5.10 and 5.11)
scenario_tasks("one_location_step")
scenario_tasks("two_location_step")
scenario_tasks("two_location_step_diff")


## In the case of "one_location_vary_si" and "one_location_vary_cv",
## we also explored a scenario where where the SI distribution
## (mean or CV) of the variant is different from that of the reference
## but in the absence of more information, is assumed to be the same
## as that of the reference. This required separate estimation
## and summary tasks.
## (See SI Sec 5.5 and 5.7)


a <- orderly_run("one_location_wrong_si",
                 parameters = list(short_run = FALSE))
orderly_commit(a)
a <- orderly_run("one_location_wrong_cv",
                 parameters = list(short_run = FALSE))
orderly_commit(a)

a <- orderly_run("summarise_one_location_wrong_si")
orderly_commit(a)
a <- orderly_run("summarise_one_location_wrong_cv")
orderly_commit(a)


## In the scenario where we assess the performance of the method for
## one location in which there is a time-varying transmission
## advantage, we also explore the impact of varying the amount of
## data we use during the estimation step. Therefore the following
## tasks have an additional parameter ("estimation_window"). We
## explore three different values of the estimation window.
## (See SI Sec 5.12)


window <- c(10, 7) #, "standard")

walk(window, function(w) {
  
  a <- orderly_run("sims_one_location_changing_advantage",
                   parameters = list(short_run = TRUE))
  orderly_commit(a)
  
  a <- orderly_run("one_location_changing_advantage",
                   parameters = list(short_run = TRUE,
                                     estimation_window = w))
  orderly_commit(a)

  a <- orderly::orderly_run("summarise_one_location_changing_advantage",
                            parameters = list(estimation_window = w))
  orderly_commit(a)

})


## Run these tasks to assess whether the new variant is more or less
## transmissible than the reference.


a <- orderly_run("classify_as_transmissible")
orderly_commit(a)
a <- orderly_run("classify_as_transmissible_step")
orderly_commit(a)
a <- orderly_run("classify_as_transmissible_changing_advantage")
orderly_commit(a)


## Run these tasks to produce figures and other results in the manuscript


# Results in main text (Fig 2) and SI Secs 5.3 - 5.9
a <- orderly_run("produce_summary_figures2")
a <- orderly_run("produce_summary_figures")

# Results in SI Secs 5.10 - 5.11
a <- orderly_run("produce_stepwise_summary_figures")

# Results in SI Secs 5.12
a <- orderly_run("produce_changing_adv_summary_figures")
