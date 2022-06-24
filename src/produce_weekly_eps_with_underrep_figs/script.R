## orderly::orderly_develop_start(use_draft = "newer")
source("R/fig_utils.R")
dir.create("figures")
whole_cnty_time <- readRDS("epsilon_qntls_whole_country.rds")[["uk_alpha_wild"]]

## Epsilon estimates with overlapping weekly windows.
## Axis 2. Proportion of variant in the window used for estimation
ovl_est <- readRDS("epsilon_qntls_over_time.rds")[["uk_alpha_wild"]]
## Epsilon estimates with overlapping weekly windows.
## and underreporting
ovl_est_underrep <- readRDS("epsilon_qntls_over_time_underep.rds")[["uk_alpha_wild"]]
## Put them together
ovl_est$date <- as.Date(ovl_est$date)
ovl_est_underrep$date <- as.Date(ovl_est_underrep$date)
ovl_est <- bind_rows(
  None = ovl_est, `50%` = ovl_est_underrep, .id = "Underreporting"
)

eng_noadj <- readRDS("england_na_not_adjusted.rds")[["alpha"]]
## Daily proportion is volatile, do weekly
eng_noadj <- slider::slide_period_dfr(eng_noadj, eng_noadj$date, "week", function(x) {
  out <- broom::tidy(apply(x[, -1], 2, sum))
  data.frame(
    date = x$date[1],
    out
  ) %>% spread(names, x)
})

out <- binconf(
  eng_noadj$wildtype + eng_noadj$alpha,
  eng_noadj$wildtype + eng_noadj$alpha + eng_noadj$unknown
)
eng_noadj <- cbind(eng_noadj, data.frame(out))

y2label <- "Proportion with variant known"
xmin <- xaxis_breaks$uk_alpha_wild
col <- palette[c("alpha")]
z <- left_join(ovl_est, eng_noadj, by = "date")
whole_cnty_time$date <- max(ovl_est$date) + 5 ## slightly to the right to avoid overlap
coeff <- 1.8 ## For everyhing except delta
z$proportion_scaled <- z$`PointEst` * coeff
z$low_scaled <- z$`Lower` * coeff
z$high_scaled <- z$`Upper` * coeff


dodge_width <- 5

p <- ggplot(z, aes(x = date)) +
  geom_point(
    aes(y = `50%`, shape = Underreporting), colour = col, size = 2,
    position = position_dodge(width = dodge_width)
  ) +
  geom_linerange(
    aes(ymin = `2.5%`, ymax = `97.5%`, linetype = Underreporting),
    size = 1.1, colour = col,
    position = position_dodge(width = dodge_width)
  ) +
  geom_hline(
    yintercept = 1, linetype = "dashed", color = "red", size = 1.1
  ) +
  geom_point(aes(y = proportion_scaled), color = "blue") +
  geom_linerange(aes(ymin = low_scaled, ymax = high_scaled), color = "blue") +
  geom_line(aes(y = proportion_scaled), color = "blue", alpha = 0.5) +
  scale_y_continuous(
    sec.axis = sec_axis(~./coeff, name = y2label, labels = mypercent),
    limits = c(0, coeff),
    breaks = seq(0, coeff, by = 0.5)
  ) +
  scale_x_date(
    breaks = xmin,
    ## date_breaks = date_breaks,
    date_labels = date_labels,
    limits = c(min(xmin), NA)
  ) +
  scale_linetype_manual(
    name = "Reporting",
    breaks = c("None", "50%"),
    labels = c("100%", "50%"),
    values = c(None = "solid", "50%" = "dashed")
  ) +
  scale_shape_manual(
    name = "Reporting",
    breaks = c("None", "50%"),
    labels = c("100%", "50%"),
    values = c(None = 19, "50%" = 1)
  ) +
  ## Add estimate for the entire country over the whole time period
  ## geom_point(
  ##   data = whole_cnty_time, aes(x = date, y = `50%`), fill = col, size = 3,
  ##   shape = 23, col = "black", stroke = 2
  ## ) +
  ## geom_linerange(
  ##   data = whole_cnty_time, aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
  ##   size = 1.1, colour = col
  ## ) +
#########
  coord_cartesian(clip = "off") +
  ylab("Effective transmission advantage") +
  xlab("") +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.line.y.right = element_line(color = "blue"),
    axis.title.y.right = element_text(color = "blue"),
    axis.text.y.right = element_text(color = "blue")
  )


save_multiple(p, glue("figures/underrep_uk_alpha_wild_2axis"))

