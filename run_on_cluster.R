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
obj <- didehpc::queue_didehpc(ctx)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
res <- obj$enqueue(orderly::orderly_bundle_run('20210713-160158-cfa36a12.zip', "cluster-runs"))
# id "3dc0bba4842f281240837c753b7db2bb"
orderly::orderly_bundle_import('cluster-runs/20210712-223938-7707dbc9.zip')
