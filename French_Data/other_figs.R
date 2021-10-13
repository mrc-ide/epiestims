library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)
theme_manuscript <- function(base_size = 14) {
  theme_minimal() %+replace%
    theme(
      text = element_text(size = base_size),
      legend.position = "top"
      ##axis.text.x = element_text(angle = 45)
    )
}
## give filename without the extension
save_multiple <- function(plot, filename) {
  ggsave(
    filename = glue("{filename}.pdf"),
    plot
  )
  ggsave(
    filename = glue("{filename}.png"),
    plot)
}


outdir <- "figs/other_figs"
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

periods <- readRDS(file = 'Rdata/periods.rds')

france <- readRDS(file = 'Rdata/I_fr.rds')
names(france) <- c("Wild", "Alpha", "Beta/Gamma")
france <- bind_rows(france, .id = "variant")

periods_fr <- periods$periods_fr$intervals
dates_fr <- france$date[periods_fr[-1]]
france <- gather(france, region, incid, -date, -variant)
france$region[france$region == "Bourgogne-Franche-Comté"] <- "Bourgogne- Franche-Comté"



p <- ggplot(france) +
  geom_line(aes(date, incid, col = variant), size = 1.1) +
  facet_wrap(
    ~region, scales = "free_y", ncol = 3,
    labeller = labeller(region = label_wrap_gen(10))
  ) +
  geom_vline(xintercept = dates_fr, linetype = "dashed") +
  scale_color_manual(
    breaks = c("Wild", "Alpha", "Beta/Gamma"),
    values = c(Wild = "#000000", Alpha = "#E69F00", `Beta/Gamma` = "#56B4E9")
  ) +
  ylab("Incidence") +
  theme_manuscript() +
  theme(axis.title.x = element_blank(), legend.title = element_blank())

save_multiple(p, glue("{outdir}/france_incidence"))

england1 <- readRDS(file = 'Rdata/I_UK1.rds')
names(england1) <- c("Wild", "Alpha")
england1 <- bind_rows(england1, .id = "variant")

periods_en <- periods$periods_UK1$intervals
dates_en <- england1$date[periods_en[-1]]
england1 <- gather(england1, region, incid, -date, -variant)


p <- ggplot(england1) +
  geom_line(aes(date, incid, col = variant), size = 1.1) +
  facet_wrap(
    ~region, scales = "free_y", ncol = 3,
    labeller = labeller(region = label_wrap_gen(10))
  ) +
  geom_vline(xintercept = dates_en, linetype = "dashed") +
  scale_color_manual(
    breaks = c("Wild", "Alpha"),
    values = c(Wild = "#000000", Alpha = "#E69F00")
  ) +
  ylab("Incidence") +
  theme_manuscript() +
  theme(axis.title.x = element_blank(), legend.title = element_blank())

save_multiple(p, glue("{outdir}/england_early_incidence"))


england2 <- readRDS(file = 'Rdata/I_UK2.rds')
names(england2) <- c("Alpha", "Delta")
england2 <- bind_rows(england2, .id = "variant")

periods_en <- periods$periods_UK2$intervals
dates_en <- england2$date[periods_en[-1]]
england2 <- gather(england2, region, incid, -date, -variant)

p <- ggplot(england2) +
  geom_line(aes(date, incid, col = variant), size = 1.1) +
  facet_wrap(
    ~region, scales = "free_y", ncol = 3,
    labeller = labeller(region = label_wrap_gen(10))
  ) +
  geom_vline(xintercept = dates_en, linetype = "dashed") +
  scale_color_manual(
    breaks = c("Alpha", "Delta"),
    values = c(Delta = "#009E73", Alpha = "#E69F00")
  ) +
  ylab("Incidence") +
  theme_manuscript() +
  theme(axis.title.x = element_blank(), legend.title = element_blank())

save_multiple(p, glue("{outdir}/england_late_incidence"))

