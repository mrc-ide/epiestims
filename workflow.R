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
