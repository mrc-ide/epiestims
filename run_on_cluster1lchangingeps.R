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
# [ init:id   ]  20c4b942489c74e3fa87c93aa92c38d1
# [ init:db   ]  rds
# [ init:path ]  context
# [ save:id   ]  e5a2e2d54b045195446a914c826dcb58
# [ save:name ]  humourful_allosaurus
# ctx <- context_read('a411c9a04e2543646549240c451d83d7', root)


a <- orderly::orderly_run("sims_one_location_changing_advantage",
                          parameters = list(short_run = FALSE))
orderly::orderly_commit(a)

path_bundles <- "Q:/cluster"

bundle1 <- orderly::orderly_bundle_pack(path_bundles,
                                        "one_location_changing_advantage",
                                        parameters = list(short_run = FALSE,
                                                          estimation_window = 10))

bundle2 <- orderly::orderly_bundle_pack(path_bundles,
                                        "one_location_changing_advantage",
                                        parameters = list(short_run = FALSE,
                                                          estimation_window = 7))

bundle3 <- orderly::orderly_bundle_pack(path_bundles,
                                        "one_location_changing_advantage",
                                        parameters = list(short_run = FALSE,
                                                          estimation_window = 1))

bundle4 <- orderly::orderly_bundle_pack(path_bundles,
                                        "one_location_changing_advantage",
                                        parameters = list(short_run = FALSE,
                                                          estimation_window = "standard"))

bundle1$path # bundle path name to include in orderly_bundle_run
bundle2$path
bundle3$path
bundle4$path
# "Q:\\cluster\\20220421-151815-d44aa9ea.zip
# Full run (all tmax values - files above are a test of two tmax)
# "Q:\\cluster\\20220421-163105-d5cdead9.zip"
# Re-run with corrected t_min values
# "Q:\\cluster\\20220509-103711-94794719.zip"
# Re-run with a 7 day time window
# "Q:\\cluster\\20220509-163536-7767a098.zip"
# Re-run with a 1 day time window
# "Q:\\cluster\\20220510-172110-1882d7d2.zip"

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim')

# orderly_bundle_run arguments
# 1) where bundle_pack was saved (find with bundle$path)
# 2) where bundle_run outputs should go (seem to need to specify full file path)
res1 <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20220610-160549-aa6c896d.zip",
                                                "Z:\\jwardle\\epiestims\\cluster-runs"))

res2 <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20220610-160556-dbc63195.zip",
                                                "Z:\\jwardle\\epiestims\\cluster-runs"))

res3 <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20220610-160602-42a25409.zip",
                                                "Z:\\jwardle\\epiestims\\cluster-runs"))

res4 <- obj$enqueue(orderly::orderly_bundle_run("Q:\\cluster\\20220610-160607-faf12c2e.zip",
                                                "Z:\\jwardle\\epiestims\\cluster-runs"))

# task id: "62f5492dc4abd0def05a2cc780812396"


orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20220610-160549-aa6c896d.zip")

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20220610-160556-dbc63195.zip")

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20220610-160602-42a25409.zip")

orderly::orderly_bundle_import("Z:\\jwardle\\epiestims\\cluster-runs\\20220610-160607-faf12c2e.zip")
