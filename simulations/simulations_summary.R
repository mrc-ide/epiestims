source("simulations/simulation_functions.R")
##source("simulations/simulations_si.R")
library(glue)
library(dplyr)
library(purrr)
library(ggplot2)

sim_params <- readRDS("results/vary_si_sim_params.rds")
prefix <- "vary_si"

fit_files <- glue(
  "results/{prefix}_{seq_len(nrow(sim_params))}.rds"
)
## we have some missing results due to cluster
## issues
available <- file.exists(fit_files)
## fit_files <- fit_files[available]
## sim_params <- sim_params[available, ]

incid_at_tmax <- imap_dfr(
  available, function(avl, index) {
    ## Can't read all the files at the same time
    ## memory exhausted
    ## read one at a time
    if (! avl) return(NULL)
    infile <- glue("results/{prefix}_{index}.rds")
    fit <- readRDS(infile)
    message("Reading ", infile)
    params <- sim_params[index, ]
    imap_dfr(
      fit$incid, function(sim, sim_index) {
        map_dfr(
        tmax_all, function(tmax) {
          out <- data.frame(
            ref_incid = sum(sim[1:tmax, 1, 1]),
            var_incid = sum(sim[1:tmax, 1, 2])
          )
          cbind(params, out)
        }, .id = "tmax")}, .id = "sim"
    )
  }
)


## Structure of each fit$res object is:
## List of 5 corresponding to the 5 tmax values
## seq(20, 60, 10)
## For each tmax value, we have a list of 100
## simulations, each of which has epsilon and R

eps_summary <- imap_dfr(
  available, function(avl, index) {
    ## Just read them again, it doesm't matter!
    ## otherwise you have  will have to make one
    ## massive loop doing everything
    if (! avl) return(NULL)
    infile <- glue("results/{prefix}_{index}.rds")
    fit <- readRDS(infile)
    message("Reading ", infile)
    params <- sim_params[index, ]
    if (is.null(fit$res)) {
      message("Result NULL at ", index)
      return(NULL)
    }
    names(fit$res) <- tmax_all
    map_dfr(
      fit$res, function(eps_tmax) {
        imap_dfr(
          eps_tmax, function(eps_sim, sim) {
            out <- summarise_epsilon(eps_sim)
            cbind(params, out)
          }, .id = "sim")
      }, .id = "tmax"
    )
  }
)




eps_err_summary <- imap_dfr(
  available, function(avl, index) {
    ##  Read them again
    if (! avl) return(NULL)
    infile <- glue("results/{prefix}_{index}.rds")
    fit <- readRDS(infile)
    message("Reading ", infile)
    params <- sim_params[index, ]
    if (is.null(fit$res)) {
      message("Result NULL at ", index)
      return(NULL)
    }
    names(fit$res) <- tmax_all
    map_dfr(
      fit$res, function(eps_tmax) {
        imap_dfr(
          eps_tmax, function(eps_sim, sim) {
            out <- summarise_epsilon_error(
              eps_sim, params$epsilon
            )
            cbind(params, out)
          }, .id = "sim")
      }, .id = "tmax"
    )
  }
)



saveRDS(
  incid_at_tmax, glue("results/incid_at_tmax_{prefix}.rds")
)

saveRDS(
  eps_err_summary, glue("results/eps_err_summary_{prefix}.rds")
)

saveRDS(
  eps_summary, glue("results/eps_summary_{prefix}.rds")
)






