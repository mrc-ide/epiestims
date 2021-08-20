library(orderly)
library(context)

options(
  # didehpc.cluster = "fi--didemrchnb",
  didehpc.cluster = "fi--dideclusthn",
  didehpc.username = "jw2519",
  didehpc.home = "Q:/")

config <- didehpc::didehpc_config(cores = 8, parallel = TRUE)

packages <- c("dplyr", "furrr", "glue", "incidence",
              "orderly", "projections", "purrr", "utils")

root <- 'context'
ctx <-context_save(root, packages = packages)
# [ init:id   ]  340a337ed2403fbbddb633d2200a2f4b
# [ init:db   ]  rds
# [ init:path ]  context
# [ save:id   ]  81a4a89d9647e2521c808e4ab78b836a
# [ save:name ]  useful_jumpingbean
# ctx <- context_read('81a4a89d9647e2521c808e4ab78b836a', root)


a <- orderly::orderly_run("sims_one_location_step",
                          parameters = list(short_run = TRUE))
b <- orderly::orderly_run("sims_two_location_step",
                          parameters = list(short_run = FALSE))
orderly::orderly_commit(b)

path_bundles <- "Q:/cluster"

bundle <- orderly::orderly_bundle_pack(path_bundles, "two_location_step",
                             parameters = list(short_run = FALSE))

bundle$path # bundle path name to include in orderly_bundle_run

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')

# orderly_bundle_run arguments
# 1) where bundle_pack was saved (find with bundle$path)
# 2) where bundle_run outputs should go (seem to need to specify full file path)
res <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20210820-152308-fe544983.zip",
                                               "Z:\\jwardle\\epiestims\\cluster-runs"))
# task id: bdb4fe069f5f908f5565af1c5b162079

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20210813-175522-5ed36934.zip")
