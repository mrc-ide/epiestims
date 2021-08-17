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
# [ init:id   ]  b02a8a928368b89edd4c4aecc42a6f58
# [ init:db   ]  rds
# [ init:path ]  context
# [ save:id   ]  542c1033c55c493e03cf574fba298f71
# [ save:name ]  churnable_nightheron
# ctx <- context_read('542c1033c55c493e03cf574fba298f71', root)


a <- orderly::orderly_run("sims_one_location_step",
                          parameters = list(short_run = TRUE))
b <- orderly::orderly_run("sims_one_location_step",
                          parameters = list(short_run = FALSE))
orderly::orderly_commit(b)

path_bundles <- "Q:/cluster"

bundle <- orderly::orderly_bundle_pack(path_bundles, "one_location_step",
                             parameters = list(short_run = FALSE))

bundle$path # bundle path name to include in orderly_bundle_run

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')

# orderly_bundle_run arguments
# 1) where bundle_pack was saved (find with bundle$path)
# 2) where bundle_run outputs should go (seem to need to specify full file path)
res <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20210813-175522-5ed36934.zip",
                                               "Z:\\jwardle\\epiestims\\cluster-runs"))
# task id: bdb4fe069f5f908f5565af1c5b162079

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20210813-175522-5ed36934.zip")
