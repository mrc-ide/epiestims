library(context)
root <- "correct_si"
packages <- c("dplyr", "here","purrr", "epitrix", "glue", "incidence", 
              "projections")
# src <- conan::conan_sources("mrc-ide/EpiEstim@multiv")
# options(didehpc.cluster = 'fi--didemrchnb')

source_files <- c(
  "global.R", "simulations/simulation_functions.R",
  "simulations/simulations_si.R",
  "simulations/cluster_functions.R"
)

ctx <-context_save(root, packages = packages, sources = source_files)
# dissimilar_teledu
# ctx <- context_read('e450caecbc8ae70f5afd841b6c616104', root)

obj <- didehpc::queue_didehpc(ctx)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@draw_eps_debug')
vary_si <- obj$enqueue_bulk(sim_params, manager)
# 27.06.2021 'dispersible_noctule'
# Mis-specification of SI task submitted 10th June 'thalassophilic_chipmunk'
## Test locally,
## obj <- queuer:::queue_local$new(ctx)
## tb <- obj$enqueue_bulk(sim_params[c(1, 11), ], manager)
## obj$run_next
# 28.05.2021 'winsome_alaskanmalamute'
vary_si <- obj$task_bundle_get('dispersible_noctule')
idx <- which(vary_si$status() == 'COMPLETE')
iwalk(
  idx[105:length(idx)],
  function(index, task) {
    t <- obj$task_get(task)
    res <- t$result()
    message('saving ', glue('results/vary_si_{index}.rds'))
    saveRDS(res, glue('results/vary_si_{index}.rds'))
  }
  )

vary_si <- obj$enqueue_bulk(sim_params, manager)
# 'antinationalistic_walrus'