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
ctx <- context_save(root, packages = packages)
# [ init:id   ]  2f08608760e2b524f962997a01ff6e13
# [ init:db   ]  rds
# [ init:path ]  context
# [ save:id   ]  4604712b0d5bae09f2a9421e7fb191b4
# [ save:name ]  unisex_vervet
# ctx <- context_read('4604712b0d5bae09f2a9421e7fb191b4', root)


a <- orderly::orderly_run("sims_one_location_changing_advantage",
                          parameters = list(short_run = FALSE))
orderly::orderly_commit(a)

path_bundles <- "Q:/cluster"

bundle1 <- orderly::orderly_bundle_pack(path_bundles,
                                        "one_location_changing_advantage",
                                        parameters = list(short_run = FALSE))

bundle1$path # bundle path name to include in orderly_bundle_run
# "Q:\\cluster\\20220407-122655-7821ec22.zip"

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim')

# orderly_bundle_run arguments
# 1) where bundle_pack was saved (find with bundle$path)
# 2) where bundle_run outputs should go (seem to need to specify full file path)
res1 <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20220407-122655-7821ec22.zip",
                                                "Z:\\jwardle\\epiestims\\cluster-runs"))

# task id: bdb4fe069f5f908f5565af1c5b162079

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20220407-122655-7821ec22.zip")
