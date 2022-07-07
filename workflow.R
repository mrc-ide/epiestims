library(orderly)

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

## Apply function to scenarios featured in our analysis

scenario_tasks("one_location_vary_si")
scenario_tasks("one_location_vary_cv")
scenario_tasks("one_location_vary_offspring")
scenario_tasks("one_location_vary_underreporting")

scenario_tasks("two_locations")
scenario_tasks("two_diff_locations")

scenario_tasks("one_location_step")
scenario_tasks("two_location_step")
scenario_tasks("two_location_step_diff")


# In the scenario where we assess the performance of the method for
# one location in which there is a time-varying transmission
# advantage, we also explore the impact of varying the amount of
# data we use during the estimation step. Therefore the following
# tasks have an additional parameter ("estimation_window").

window <- c(10, 7, "standard")

map(window, function(w) {

a <- orderly_run(sims_one_location_changing_advanatge,
                 parameters = list(short_run = FALSE))
orderly_commit(a)

a <- orderly_run(one_location_changing_advantage,
                 parameters = list(short_run = FALSE,
                                   estimation_window = w))
orderly_commit(a)

a <- orderly::orderly_run("summarise_one_location_changing_advantage",
                          parameters = list(estimation_window = w))
orderly_commit(a)

})



## We can first generate the simulated data for a range of
## scenarios with the following orderly tasks.

a <- orderly_run("sims_one_location_vary_si")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_cv")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_offspring")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_underreporting")
orderly_commit(a)

a <- orderly_run("sims_one_location_step")
orderly_commit(a)

a <- orderly_run("sims_two_locations") ##NOTE: DID WE USE THE 2 LOCATION SCENARIOs??
orderly_commit(a)

a <- orderly_run("sims_two_diff_location") 
orderly_commit(a)

a <- orderly_run("sims_two_location_step")
orderly_commit(a)

a <- orderly_run("sims_two_location_step_diff")
orderly_commit(a)

a <- orderly_run("sims_one_location_changing_advantage")
orderly_commit(a)


## We then apply our method to estimate the transmission
## advantage on each set of simulated incidence data.

a <- orderly_run("sims_one_location_vary_si")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_cv")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_offspring")
orderly_commit(a)

a <- orderly_run("sims_one_location_vary_underreporting")
orderly_commit(a)

a <- orderly_run("sims_one_location_step")
orderly_commit(a)

a <- orderly_run("sims_two_locations") ##NOTE: DID WE USE THE 2 LOCATION SCENARIOs??
orderly_commit(a)

a <- orderly_run("sims_two_diff_location") 
orderly_commit(a)

a <- orderly_run("sims_two_location_step")
orderly_commit(a)

a <- orderly_run("sims_two_location_step_diff")
orderly_commit(a)

a <- orderly_run("sims_one_location_changing_advantage")
orderly_commit(a)



a <- orderly_run("summarise_one_location_vary_si")
orderly_commit(a)

a <- orderly_run("summarise_one_location_wrong_si")
orderly_commit(a)

a <- orderly_run("summarise_one_location_vary_cv")
orderly_commit(a)

a <- orderly_run("summarise_one_location_wrong_cv")
orderly_commit(a)

a <- orderly_run("summarise_one_location_vary_offspring")
orderly_commit(a)

a <- orderly_run("summarise_one_location_underreporting")
orderly_commit(a)

a <- orderly_run("classify_as_transmissible")
orderly_commit(a)



a <- orderly_run("produce_summary_figures2")
a <- orderly_run("produce_summary_figures")
## Real data

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


