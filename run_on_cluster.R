library(context)
packages <- c("dplyr", "epitrix", "glue", "incidence",
              "orderly", "projections", "purrr")
obj <- didehpc::queue_didehpc(ctx)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
res <- obj$enqueue(orderly::orderly_bundle_run(, "cluster-runs"))
