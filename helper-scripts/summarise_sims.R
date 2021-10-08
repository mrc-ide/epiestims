## This is common to all summary tasks
dir.create("figures")
zip::unzip("output.zip")
sim_params <- readRDS("param_grid.rds")
incid <- readRDS("incid.rds")


## Summarise simulated incidence, summarise
## epsilon and epsilon error.
tmax_all <- seq(10, 50, by = 10)
names(tmax_all) <- tmax_all
si_mu_ref <- 5.4
si_std_ref <- 1.5
si_distr_ref <- discr_si(0:30, mu = si_mu_ref, sigma = si_std_ref)
incid_summary <- map_depth(
  incid, 2, function(x) {
    map_dfr(
      tmax_all, function(tmax) {
        data.frame(
          ref_incid = sum(x[1:tmax, 1, 1]),
          var_incid = sum(x[1:tmax, 1, 2]),
          tmax = tmax
        )
      }
    )
  }
)

incid_summary <- map(
  incid_summary, ~ bind_rows(., .id = "sim")
)

saveRDS(incid_summary, "incid_summary.rds")

## Structure of output: list with one element
## for each row of params.
## Each element is a list of length 5, corresponding
## to the 5 tmax values we use (tmin is set to 15).
## Each element at nested level 2 is a list of
## length 100, corresponding to the number of
## simulations that we run.
## Finally within that, we have a list of length
## 2, which is the output from estimate_joint.
eps_summary <- map(
  seq_len(nrow(sim_params)),
  function(index) {
    infile <- glue("outputs/estimate_joint_{index}.rds")
    if (!file.exists(infile)) {
      warning(infile, " not present")
      return(NULL)
    }
    res <- readRDS(infile)
    message("At ", infile)
    names(res) <- tmax_all
    out <- map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon(res_sim[[1]])
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
    out <- estimate_uncertain(out)
  }
)

## Note Pierre's method is too strict.
## From Pierre's method, we get something like this:
## Obviously cases we don't want to classify as having
## high uncertainty.
## > tail(low)
##       tmax sim  2.5%   25%   50%   75% 97.5%       mu        sd   param rt_ref
## 51700   20 100 2.556 2.859 3.006 3.162 3.520 3.009315 0.2394802 epsilon    1.6
## 51706   30   6 2.744 2.896 2.982 3.053 3.245 2.981690 0.1245929 epsilon    1.6
## 51747   30  47 2.685 2.842 2.917 2.998 3.189 2.920902 0.1227537 epsilon    1.6
## 51755   30  55 2.438 2.609 2.700 2.800 2.957 2.702453 0.1380378 epsilon    1.6
## 51763   30  63 2.427 2.588 2.686 2.778 2.939 2.689514 0.1304141 epsilon    1.6
## 51768   30  68 2.648 2.812 2.905 3.000 3.166 2.908238 0.1367992 epsilon    1.6
##       true_eps si_mu_variant cri_width confidence                 true_label
## 51700        3          10.8     0.964        Low Variant more transmissible
## 51706        3          10.8     0.501        Low Variant more transmissible
## 51747        3          10.8     0.504        Low Variant more transmissible
## 51755        3          10.8     0.519        Low Variant more transmissible
## 51763        3          10.8     0.512        Low Variant more transmissible
## 51768        3          10.8     0.518        Low Variant more transmissible

eps_err_summary <- map2(
  seq_len(nrow(sim_params)),
  sim_params$epsilon,
  function(index, epsilon) {
    infile <- glue("outputs/estimate_joint_{index}.rds")
    if (!file.exists(infile)) {
      warning(infile, " not present")
      return(NULL)
    }
    res <- readRDS(infile)
    message("At ", infile)
    names(res) <- tmax_all
    map_dfr(
      res, function(res_tmax) {
        map_dfr(
          res_tmax, function(res_sim) {
            summarise_epsilon_error(res_sim[[1]], epsilon)
          }, .id = "sim"
        )
      }, .id = "tmax"
    )
  }
)
