library(orderly)
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
# [ init:id   ]  6b39f6628ea252c37ab705e9aa91861c
# [ init:db   ]  rds
# [ init:path ]  context
# [ save:id   ]  f1d1ce5a48c80487db1dfe49a7e2c8b0
# [ save:name ]  centrifugal_paintedladybutterfly
# ctx <- context_read('f1d1ce5a48c80487db1dfe49a7e2c8b0', root)


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
res <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20210817-114617-7314fda5.zip",
                                               "Z:\\jwardle\\epiestims\\cluster-runs"))
# task id: bdb4fe069f5f908f5565af1c5b162079

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20210813-175522-5ed36934.zip")

