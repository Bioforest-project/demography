graph_traj <- function(data, lab_y) {
  controls <- data |>
    group_by(Site) |>
    mutate(tmin = min(trec[type == "rec"]), tmax = max(trec[type == "rec"])) |>
    filter(type == "equ") |>
    mutate(trec = ifelse(Treatment == "Logging", tmin - 2, tmax + 2))
  data |>
    filter(type == "rec") |>
    ggplot(aes(trec, y = obs)) +
    geom_line(aes(group = Plot, col = Plot, y = pred), size = 1.2) +
    geom_point(aes(group = Plot, col = Plot)) +
    geom_line(aes(group = Plot, col = Plot)) +
    geom_boxplot(data = controls, aes(group = Treatment), width = 2) +
    xlab("Recovery time [yr]") +
    ylab(lab_y) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(~Site, scales = "free")
}
