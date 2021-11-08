library(orderly)

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



a <- orderly_run("produce_summary_figures")

## Real data

a <- orderly_run("clean_french_england_data")
orderly_commit(a)

a <- orderly_run("naive_epsilon_estimates")
orderly_commit(a)

a <- orderly_run("mv_epiestim_weekly_estimates")
orderly_commit(a)


a <- orderly_run("mv_epiestim_regional_estimates")
orderly_commit(a)

