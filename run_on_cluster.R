config <- didehpc::didehpc_config(cores = 8, parallel = TRUE)
library(context)
packages <- c("dplyr", "epitrix", "glue", "incidence",
              "orderly", "projections", "purrr")
root <- 'context'
ctx <-context_save(root, packages = packages)
# [ open:db   ]  rds
# [ save:id   ]  ed59bb10737dc761b869e66663baef1d
# [ save:name ]  unadored_lemur
ctx <- context_read('ed59bb10737dc761b869e66663baef1d', root)
# orderly::orderly_bundle_pack(".", "one_location_vary_si", parameters = list(short_run = FALSE))
obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@draw_eps_debug')
res <- obj$enqueue(orderly::orderly_bundle_run('20210716-121320-d1cfd238.zip', "cluster-runs"))
# task id: 7d851b4088f1abe608d50ebd78d865ca
