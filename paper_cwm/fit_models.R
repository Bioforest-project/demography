library(tidyverse)
library(cmdstanr)
library(bayesplot)

setwd("paper_cwm")

source("functions.R")
data <- read.csv("data/data.csv")
data_rec <- data |>
  filter(Treatment == "Logging" & cohort != "all") |>
  mutate(group = paste0(cohort, "_", flux)) |>
  mutate(sitenum = as.numeric(as.factor(Site))) |>
  mutate(plotnum = as.numeric(as.factor(paste(Site, plot_new)))) |>
  mutate(traitnum = as.numeric(as.factor(trait))) |>
  mutate(trec = Year - Year_of_harvest) |>
  mutate(t0 = t0 - Year_of_harvest) |>
  # cumulative changes
  group_by(Site, plot_new, trait, group) |>
  arrange(Year) |>
  mutate(cum_y = cumsum(value)) |>
  ungroup() |>
  # temp: remove SUAS
  filter(Site != "SUAS")

groups <- unique(data_rec$group)

## annual model
for (gr in groups) {
  data_group <- data_rec |> filter(group == gr)
  mdata_cwd <- list(
    N = nrow(data_group),
    S = max(data_group$sitenum),
    P = max(data_group$plotnum),
    K = max(data_group$traitnum),
    cum_y = data_group$cum_y,
    trec = data_group$trec,
    t0 = data_group$t0,
    site = data_group$sitenum,
    plot = data_group$plotnum,
    trait = data_group$traitnum,
    tau_range = c(5, 25),
    gamma_range = c(1, 10)
  )
  sample_model("annual_bump", mdata_cwd, paste0("annual_bump_", gr))
}

## cumsum model
for (gr in groups) {
  data_group <- data_rec |> filter(group == gr)
  mdata_cwd <- list(
    N = nrow(data_group),
    S = max(data_group$sitenum),
    P = max(data_group$plotnum),
    K = max(data_group$traitnum),
    cum_y = data_group$cum_y,
    trec = data_group$trec,
    t0 = data_group$t0,
    site = data_group$sitenum,
    plot = data_group$plotnum,
    trait = data_group$traitnum,
    tau_range = c(5, 25),
    gamma_range = c(1, 10)
  )
  sample_model("cumulative_bump", mdata_cwd, paste0("cumulative_bump_", gr))
}


## LOO bump selection
