config <- didehpc::didehpc_config(cores = 8, parallel = TRUE)
library(context)
packages <- c("dplyr", "epitrix", "furrr" ,"glue", "incidence",
              "orderly", "projections", "purrr")
root <- 'vary_si'
ctx <-context_save(root, packages = packages)
# [ open:db   ]  rds
# [ save:id   ]  c061ec7f3677d8e8765c86e53326c716
# [ save:name ]  quasiobjective_indianspinyloach
ctx <- context_read('c061ec7f3677d8e8765c86e53326c716', root)
# orderly::orderly_bundle_pack(".", "one_location_vary_si", parameters = list(short_run = FALSE))
obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
res <- obj$enqueue(orderly::orderly_bundle_run("Z:\\sbhatia\\epiestims\\20210813-183738-87e2789a.zip", "Z:\\sbhatia\\epiestims\\cluster-runs"))
# task id: 3623f0cea161fbee7fa60bab1a006bfe
