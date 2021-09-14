dir.create("figures")

## Function for saving
save_multiple <- function(plot, filename) {
  ggsave(
    filename = glue("{filename}.pdf"),
    plot
  )
  ggsave(
    filename = glue("{filename}.png"),
    plot)
}

## two_location_step scenario
## generate dtaframe of Rt over time

two_loc_step <- data.frame("time" = rep(seq(1,60,1), 4),
                           "rt" = c(rep(1.4, 19), rep(1.1, 41),
                                    rep(1.4, 39), rep(1.1, 21),
                                    rep(1.6, 19), rep(1.2, 41),
                                    rep(1.6, 39), rep(1.2, 21)
                           ),
                           "location" = rep(c(rep("A", 60),
                                              rep("B", 60)),
                                            2),
                           "scenario" = c(rep(1, 120), rep(2, 120))
                            
)

two_loc_step$location <- as.factor(two_loc_step$location)
two_loc_step$scenario <- as.factor(two_loc_step$scenario)

## plot dataframe

p1 <-
  ggplot(two_loc_step, aes(x = time, y = rt)) +
  geom_line(aes(colour = location, linetype = scenario), size = 1) +
  labs(x = "Time (days)", y = "Rt", colour = "Location", linetype = "Scenario") +
  ylim(c(0, 2)) +
  theme_minimal() +
  theme(
    text = element_text(size = base_size),
    legend.position = "top"
  )

save_multiple(p1, "figures/rt_two_location_step")


## two_location_step_diff scenario
## generate dtaframe of Rt over time

two_loc_step_diff <- data.frame("time" = rep(seq(1,60,1), 2),
                           "rt" = c(rep(1.4, 39), rep(1.1, 21),
                                    rep(1.6, 19), rep(1.2, 41)
                           ),
                           "location" = c(rep("A", 60),
                                              rep("B", 60))
                           
)

two_loc_step_diff$location <- as.factor(two_loc_step_diff$location)

## plot dataframe

p2 <-
ggplot(two_loc_step_diff, aes(x = time, y = rt)) +
  geom_line(aes(colour = location), size = 1) +
  labs(x = "Time (days)", y = "Rt", colour = "Location") +
  ylim(c(0, 2)) +
  theme_minimal() +
  theme(
    text = element_text(size = base_size),
    legend.position = "top"
  )

save_multiple(p2, "figures/rt_two_location_step_diff")


