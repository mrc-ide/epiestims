library(context)
root <- "context2"
packages <- c("dplyr", "here","purrr", "epitrix", "glue", "incidence", 
              "projections")
# src <- conan::conan_sources("mrc-ide/EpiEstim@multiv")


source_files <- c(
  "global.R", "simulations/simulation_functions.R",
  "simulations/simulations_si.R",
  "simulations/cluster_functions.R"
)

ctx <-context_save(
  root, packages = packages, sources = source_files
)



obj <- didehpc::queue_didehpc(ctx)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
vary_si <- obj$enqueue_bulk(sim_params, manager)
## Test locally,
## obj <- queuer:::queue_local$new(ctx)
## tb <- obj$enqueue_bulk(sim_params[c(1, 11), ], manager)
## obj$run_next
# 28.05.2021 'winsome_alaskanmalamute'
iwalk(
  vary_si$results(), 
  function(res, index) saveRDS(res, glue('results/vary_si_{index}.rds')))

vary_si <- obj$enqueue_bulk(sim_params, manager)
# 'antinationalistic_walrus'