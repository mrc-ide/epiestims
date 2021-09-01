library(orderly)
library(context)

config <- didehpc::didehpc_config(cores = 8, parallel = TRUE)

packages <- c("dplyr", "epitrix", "furrr", "glue", "incidence",
              "projections", "purrr")

root <- 'underreport'
# ctx <-context_save(root, packages = packages)
#[ init:id   ]  4b34ac5ee19b12d136267421b8027cd0
#[ init:db   ]  rds
#[ init:path ]  underreport
#[ save:id   ]  158bf86a28432169c9de38d6a2a2dd48
#[ save:name ]  batty_hogget
# Line number 13 is Needed next time you login to the cluster
 ctx <- context_read('158bf86a28432169c9de38d6a2a2dd48', root)

# pack your task so that it can be run on the cluster:
orderly::orderly_bundle_pack(".", "sims_one_location_underreporting", 
                             parameters = list(short_run = FALSE))
#$id
#[1] "20210824-171904-1c0b630e"

#$path
#[1] "Z:\\rnash\\epiestims\\20210824-171904-1c0b630e.zip"

orderly::orderly_bundle_pack(".", "one_location_underreporting", 
                             parameters = list(short_run = FALSE))

#$id
#[1] "20210825-084703-0b03952d"

#$path
#[1] "Z:\\rnash\\epiestims\\20210825-084703-0b03952d.zip"

obj <- didehpc::queue_didehpc(ctx, config = config)
obj$install_packages('abind')
obj$install_packages('mrc-ide/EpiEstim@multiv')
obj$install_packages("vimc/orderly@fix-zip-list")

# Now run your task on the server
# 1) first argument is the $path from orderly_bundle_pack
# 2) second argument is the dir where you want orderly to write outputs
t1 <- obj$enqueue(orderly::orderly_bundle_run("Z:\\rnash\\epiestims\\20210824-171904-1c0b630e.zip", 
                                              "Z:\\rnash\\epiestims\\cluster-runs"))
# t1$id "64376932b71466e78d4e933bf77b591c"
# check t1$status()
# when t1$status() is complete, look for t1$result() and note id and path

#$id
#[1] "20210824-171904-1c0b630e"
#$path
#[1] "Z:\\rnash\\epiestims\\cluster-runs\\20210824-171904-1c0b630e.zip"


t2 <- obj$enqueue(orderly::orderly_bundle_run("Z:\\rnash\\epiestims\\20210825-084703-0b03952d.zip", 
                                              "Z:\\rnash\\epiestims\\cluster-runs"))

# t2$id "14c2464a74bac2970a0694ac5fe9021b"
# check t2$status()
# when t2$status() is 'complete', look for t2$result() and note id and path
#$id
#[1] "20210825-084703-0b03952d"
#$path
#[1] "Z:\\rnash\\epiestims\\cluster-runs\\20210825-084703-0b03952d.zip"

# now run orderly::orderly_bundle_import(<path from t1$result()>)
orderly::orderly_bundle_import("Z:\\rnash\\epiestims\\cluster-runs\\20210824-171904-1c0b630e.zip")

orderly::orderly_bundle_import("Z:\\rnash\\epiestims\\cluster-runs\\20210825-084703-0b03952d.zip")


