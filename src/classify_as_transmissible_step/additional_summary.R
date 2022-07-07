## orderly::orderly_develop_start()
source("R/utils.R")
round_to <- 3

#################################################
#################################################
####### ONE LOCATION STEPWISE #######
#################################################
#################################################

one_loc_step_eps_summary <- readRDS("one_loc_step_eps_summary_df.rds")
# one_loc_step_eps_summary <- one_loc_step_eps_summary[one_loc_step_eps_summary$confidence != "Low", ]
one_loc_step_eps_summary <- mutate_at(
  one_loc_step_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
one_loc_step_eps_summary$true_eps <- round(
  one_loc_step_eps_summary$true_eps, round_to
)
## What should be the correct classification,
## based on true advantage.
one_loc_step_eps_summary$true_label <- true_class(one_loc_step_eps_summary)
classified <- classify_epsilon(one_loc_step_eps_summary)
x <- select(classified, tmax, sim, rt_ref:est_class)
tall <- summary_tmax_eps_stepwise(x)
saveRDS(tall, "one_loc_step_classified.rds")

## Summary across rt change
classified$rt_change <- paste(classified$rt_ref, classified$rt_post_step,
                              sep = " -> ")
y <- split(classified, classified$rt_change) %>%
  map_dfr(summary_other, .id = "Change in Rt")

cat(
  stargazer(y, summary = FALSE),
  file = "one_loc_step_by_rtchange_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "one_loc_step_by_tmax_classification.tex"
)


#################################################
#################################################
####### TWO LOCATION STEPWISE #######
#################################################
#################################################
two_loc_step_eps_summary <- readRDS("two_loc_step_eps_summary_df.rds")
# one_loc_step_eps_summary <- one_loc_step_eps_summary[one_loc_step_eps_summary$confidence != "Low", ]
two_loc_step_eps_summary <- mutate_at(
  two_loc_step_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
two_loc_step_eps_summary$true_eps <- round(
  two_loc_step_eps_summary$true_eps, round_to
)
## What should be the correct classification,
## based on true advantage.
two_loc_step_eps_summary$true_label <- true_class(two_loc_step_eps_summary)
classified <- classify_epsilon(two_loc_step_eps_summary)
x <- select(classified, tmax, sim, rt_ref_l1:est_class)
tall <- summary_tmax_eps_stepwise(x)
saveRDS(tall, "two_loc_step_classified.rds")

## Summary across rt change
classified$rt_change <- paste(classified$rt_ref_l1, classified$rt_post_step_l1,
                              sep = " -> ")
y <- split(classified, classified$rt_change) %>%
  map_dfr(summary_other, .id = "Change in Rt")

cat(
  stargazer(y, summary = FALSE),
  file = "two_loc_step_by_rtchange_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "two_loc_step_by_tmax_classification.tex"
)

#################################################
#################################################
####### TWO LOCATION STEPWISE (DIFFERENT Rt) #######
#################################################
#################################################
two_loc_step_diff_eps_summary <- readRDS("two_loc_step_diff_eps_summary_df.rds")
# one_loc_step_eps_summary <- one_loc_step_eps_summary[one_loc_step_eps_summary$confidence != "Low", ]
two_loc_step_diff_eps_summary <- mutate_at(
  two_loc_step_diff_eps_summary, vars(`2.5%`:`97.5%`),
  round, round_to
)
two_loc_step_diff_eps_summary$true_eps <- round(
  two_loc_step_diff_eps_summary$true_eps, round_to
)
## What should be the correct classification,
## based on true advantage.
two_loc_step_diff_eps_summary$true_label <- true_class(two_loc_step_diff_eps_summary)
classified <- classify_epsilon(two_loc_step_diff_eps_summary)
x <- select(classified, tmax, sim, rt_ref_l1:est_class)
tall <- summary_tmax_eps_stepwise(x)
saveRDS(tall, "two_loc_step_diff_classified.rds")

## Summary across rt change
loc1_change <- paste(classified$rt_ref_l1,
                     classified$rt_post_step_l1,
                     sep = " -> ")
loc2_change <- paste(classified$rt_ref_l2,
                     classified$rt_post_step_l2,
                     sep = " -> ")


classified$rt_change <- paste("Location 1: ", loc1_change,
                              ", Location 2: ", loc2_change)
                              
y <- split(classified, classified$rt_change) %>%
  map_dfr(summary_other, .id = "Change in Rt")

cat(
  stargazer(y, summary = FALSE),
  file = "two_loc_step_diff_by_rtchange_classification.tex"
)

## Summary across tmax
y <- split(classified, classified$tmax) %>%
  map_dfr(summary_other, .id = "tmax")

cat(
  stargazer(y, summary = FALSE),
  file = "two_loc_step_diff_by_tmax_classification.tex"
)