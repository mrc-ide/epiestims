library(context)
root <- "context"
packages <- c("dplyr","purrr", "epitrix", "glue")
src <- conan::conan_sources("mrc-ide/EpiEstim@multiv")
ctx <- context_save(root, packages = packages, package_sources = src)

source_files <- c(
  "global.R", "simulations/simulation_functions.R",
  "simulations/simulations_si.R",
  "simulations/cluster_functions.R"
)

ctx <-context_save(
  root, packages = packages, sources = source_files
)



obj <- didehpc::queue_didehpc(ctx)
vary_si <- obj$enqueue_bulk(sim_params, manager)
## Test locally,
## obj <- queuer:::queue_local$new(ctx)
## tb <- obj$enqueue_bulk(sim_params[c(1, 11), ], manager)
## obj$run_next

