library(context)
packages <- c("dplyr", "epitrix", "glue", "incidence",
              "orderly", "projections", "purrr")
root <- 'context'
ctx <-context_save(root, packages = packages)
# [ open:db   ]  rds
# [ save:id   ]  ed59bb10737dc761b869e66663baef1d
# [ save:name ]  unadored_lemur
obj <- didehpc::queue_didehpc(ctx)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
res <- obj$enqueue(orderly::orderly_bundle_run('20210712-223938-7707dbc9.zip', "cluster-runs"))
# id "4605da89c0e0ea7451c7461bea32c24c"
orderly::orderly_bundle_import('cluster-runs/20210712-223938-7707dbc9.zip')
