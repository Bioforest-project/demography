library(ggplot2)
library(dplyr)

set.seed(123)

# ---- Canvas: pixel coordinates matching the original (1024 x 576) -----------
H <- 576
fy <- function(y) H - y # flip y so we can think in "image" coordinates
col_txt <- "grey30"

# ---- Simulate the community --------------------------------------------------
n_t <- 34
comm <- tibble(
  id    = seq_len(n_t),
  x     = runif_in(n_t, box_t$x),
  y     = fy(runif_in(n_t, box_t$y)),
  trait = rnorm(n_t, 1, 0.2),
  dbh   = rlnorm(n_t, 1, 1)
)

dead <- comm %>%
  arrange(trait) %>%
  slice_tail(n = 8)
survivors <- comm %>%
  filter(!id %in% dead$id) %>%
  mutate(dbh_new = dbh + rlnorm(n_t - nrow(dead), 0.2 * trait, 0.01))
recruits <- tibble(
  x = runif_in(10, box_t$x) + dx,
  y = fy(runif_in(10, box_t$y)),
  trait = rnorm(10, 1.2, 0.5),
  dbh = rnorm(10, min(comm$dbh), 0.1 * min(comm$dbh))
)

# ---- Assemble all points to draw --------------------------------------------
cloud <- function(df, xr, yr, trait = df$trait) {
  df %>% mutate(
    x = runif(n(), xr[1], xr[2]),
    y = fy(runif(n(), yr[1], yr[2])),
    trait = trait
  )
}

pts <- bind_rows(
  # community at t
  comm %>% select(x, y, trait, dbh),
  # community at t+1: survivors (grown) + recruits
  survivors %>% transmute(dbh = dbh_new),
  recruits,
  # process clouds
  cloud(recruits, c(298, 420), c(28, 135)),
  cloud(dead, c(298, 425), c(412, 525)),
  cloud(survivors, c(535, 655), c(75, 200))
)

# ---- Loop arrow for "Growth" -------------------------------------------------
theta <- seq(pi * 0.95, -pi * 0.55, length.out = 100) # clockwise sweep
loop <- tibble(
  x = 585 + 30 * cos(theta),
  y = fy(275) + 27 * sin(theta)
)

# ---- Plot -------------------------------------------------------------------
ggplot() +
  # boxes
  annotate("rect",
    xmin = box_t$x[1], xmax = box_t$x[2],
    ymin = fy(box_t$y[2]), ymax = fy(box_t$y[1]),
    fill = NA, colour = col_txt, linewidth = 0.4
  ) +
  annotate("rect",
    xmin = box_t$x[1] + dx, xmax = box_t$x[2] + dx,
    ymin = fy(box_t$y[2]), ymax = fy(box_t$y[1]),
    fill = NA, colour = col_txt, linewidth = 0.4
  ) +

  # individuals
  geom_point(
    data = pts, aes(x, y, size = dbh, colour = trait),
    show.legend = FALSE
  ) +
  scale_colour_viridis_c(option = "plasma") +

  # main arrow (t -> t+1)
  annotate("segment",
    x = 262, xend = 745, y = fy(309), yend = fy(309),
    linewidth = 2.2, colour = "black",
    arrow = arrow(length = unit(0.22, "inches"), type = "closed")
  ) +

  # recruitment and mortality arrows
  annotate("segment",
    x = 346, xend = 366, y = fy(166), yend = fy(288),
    linewidth = 0.9, colour = col_txt,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed")
  ) +
  annotate("segment",
    x = 360, xend = 380, y = fy(310), yend = fy(430),
    linewidth = 0.9, colour = col_txt,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed")
  ) +

  # growth loop
  geom_path(
    data = loop, aes(x, y), colour = col_txt, linewidth = 1.2,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed")
  ) +

  # labels
  annotate("text",
    x = 146, y = fy(168), label = "Community at t",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 857, y = fy(168), label = "Community at t+1",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 120, y = fy(465), label = "CWT = 1",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 863, y = fy(447), label = "CWT = 1.6",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 160, y = fy(38), label = "Recruitment",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 180, y = fy(70), label = "\u0394CWT = +0.2",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 597, y = fy(218), label = "Growth",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 574, y = fy(332), label = "\u0394CWT = +0.1",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 487, y = fy(490), label = "Mortality",
    size = 5.6, colour = col_txt
  ) +
  annotate("text",
    x = 510, y = fy(525), label = "\u0394CWT = +0.3",
    size = 5.6, colour = col_txt
  ) +
  coord_fixed(xlim = c(0, 1024), ylim = c(0, H), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("paper_cwm/figures/cwt_illustration.png",
  width = 10.24, height = 5.76, dpi = 100,
  bg = "white"
)
