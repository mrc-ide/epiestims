## df is a summary data.frame with columns
## 2.5% and 97.5%
classify_epsilon <- function(df) {
  df$est_class <- case_when(
    df$`2.5%` > 1 ~ "Variant more transmissible",
    df$`97.5%` < 1 ~ "Variant less transmissible",
    TRUE ~ "Unclear"
  )
  df
}

true_class <- function(df) {
  case_when(
    ## In our case this is when eps is actually 1
    df$true_eps <= 1 ~ "No transmission advantage",
    ## All other scenarios are "More transmissible"
    ## so the CrI should comfortably exclude 1.
    TRUE ~ "Variant more transmissible"
  )
}

## Summarise classification performance as a
## function of tmax and true advantage.
summary_tmax_eps <- function(x) {
  x <- split(x, x$confidence) %>%
    map_dfr(
      function(y) {
        tabyl(y, true_eps, est_class, tmax) %>%
          adorn_percentages("row") %>%
          bind_rows(.id = "tmax")
      },.id = "confidence"
    )
  
  gather(
    x, classification, val, `Unclear`:`Variant more transmissible`
  )
}

## Summarise classification performance as a
## function of tmax and true advantage (when x has no confidence variable).
summary_tmax_eps_stepwise <- function(x) {
  x <- tabyl(x, true_eps, est_class, tmax) %>%
    adorn_percentages("row") %>%
    bind_rows(.id = "tmax")
  
  gather(
    x, classification, val, `Unclear`:`Variant more transmissible`
  )
}


## More coarse summary
## x is a subset of the classified data.frame
## holding a set of variables constant.
summary_other <- function(x) {
  tabyl(x, true_label, est_class) %>%
    adorn_percentages("row")
}