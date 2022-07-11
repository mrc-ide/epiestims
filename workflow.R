library(orderly)
library(purrr)


###########################
## Estimating the effective transmission advantage
## of VOCs over the co-circulating variants using incidence
## from England and France.
## Note that this script describes a interconnected series of 'tasks' or
## analyses. This means that tasks must be run in series.
###########################

## Analysis, results and figures for main text (Fig 1) and
## SI Secs 2 & 4
## This task cleans and formats the incidence data from France and England
## for use in downstream analysis.
a <- orderly_run("clean_french_england_data")
## orderly_commit moves the task from 'draft' folder to 'archive' folder.
## Refer to orderly documentation on CRAN for more details.
## Please ensure that you run orderly_commit after each task
## so that the outputs are available for use by downstream tasks.
orderly_commit(a)

## This task will produce the non-parameteric estimates of the
## effective transmission advantage.
a <- orderly_run("naive_epsilon_estimates")
orderly_commit(a)

## This task estimates the transmission advantage of the VOcs
## (a) using only the latest 7 days of data, and (b) using the entire time
## series up to the point of estimation.
## Output RDS files containing estimatines from weekly windows are prefixed with 'nonoverlapping'.
## The window over which
## epsilon is estimated can be modfied by changing the variable
## "window" (Line 105) of the script weekly_estimates.R
a <- orderly_run("mv_epiestim_weekly_estimates")
orderly_commit(a)
## Same as above but using 50\% of cases of Alpha
## and Wildtype. This is only done for Alpha and in England.
a <- orderly_run("mv_epiestim_weekly_estimates_underrep")
orderly_commit(a)

## This task estimates the transmission advantage of VOC
## independently for each region, assuming first that there is no
## temporal variability. That is estimates are based on all the data available
## for the region up to the point of estimation. Next, we assume that
## that the effective transmission advantage also varies over time and
## use only the latest 7 days of data for estimation.
a <- orderly_run("mv_epiestim_regional_estimates")
orderly_commit(a)

## The following tasks produce summary figures and use the
## outputs from the tasks above.

## Panel A in Fig 1 and SI Figs S6, S8, and S10.
## Figs S1 and S2
orderly_run("src/produce_incid_figures")

## Panel B in Fig 1 and SI Figs S6, S8, and S10.
orderly_run("src/produce_naive_estim_figures")

## Figure S4
orderly_run("src/produce_weekly_eps_with_underrep_figs")

## Panel C in Fig 1 and SI Figs S6, S8, and S10.
a <- orderly_run("src/produce_realdata_figures/")

## Panel D in Fig 1 and SI Figs S6, S8, and S10.
## Figs S5, S7, S9, and S11.
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

## NOTE: due to estimating the transmission advantage across 100
## epidemic trajectories for each scenario and across a range of
## parameter values, the estimation_task step can take a long
## time to run. We reduced the run time of our analysis by using
## parallel computing on a cluster.
## You can edit the default parameters of the simulation, including the number
## of simulated data set for each parameter combination here
## https://github.com/mrc-ide/epiestims/blob/228ba44f0a12fb4e16346ddfaf1cee4dc9ab24b2/helper-scripts/sim_utils.R#L84
## Alternatively, you can edit the tasks that simulate data set.
## These begin with the word 'sims' e.g., sims_one_location_vary_cv


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
