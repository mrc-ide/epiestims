library(orderly)
library(context)

config <- didehpc::didehpc_config(cores = 8, parallel = TRUE)

packages <- c("dplyr", "epitrix", "furrr", "glue", "incidence",
              "projections", "purrr")

root <- 'underreport'
ctx <-context_save(root, packages = packages)
#[ init:id   ]  4b34ac5ee19b12d136267421b8027cd0
#[ init:db   ]  rds
#[ init:path ]  underreport
#[ save:id   ]  158bf86a28432169c9de38d6a2a2dd48
#[ save:name ]  batty_hogget
# Line number 13 is Needed next time you login to the cluster
# ctx <- context_read('158bf86a28432169c9de38d6a2a2dd48', root)

# pack your task so that it can be run on the cluster:
orderly::orderly_bundle_pack(".", "sims_one_location_underreporting", 
                             parameters = list(short_run = FALSE))
#$id
#[1] "20210823-165536-ac28064f"

#$path
#[1] "Z:\\rnash\\epiestims\\20210823-165536-ac28064f.zip"

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@return_diag')
obj$install_packages("vimc/orderly@fix-zip-list")

# Now run your task on the server
# 1) first argument is the $path from orderly_bundle_pack
# 2) second argument is the dir where you want orderly to write outputs
t1 <- obj$enqueue(orderly::orderly_bundle_run("Z:\\rnash\\epiestims\\20210823-165536-ac28064f.zip", 
                                              "Z:\\rnash\\epiestims\\cluster-runs"))
# t1$id b1d05d9ea2473d57a687ac1a56f6a88f
# check t1$status()
# when t1$status() is complete, look for t1$result()

#$id
#[1] "20210823-165536-ac28064f"
#$path
#[1] "Z:\\rnash\\epiestims\\cluster-runs\\20210823-165536-ac28064f.zip"

# now run orderly::orderly_bundle_import(<path from t1$result()>)
orderly::orderly_bundle_import("Z:\\rnash\\epiestims\\cluster-runs\\20210823-165536-ac28064f.zip")


