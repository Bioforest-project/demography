library(ggplot2)
library(tidyr)

visualisation <- function(y, p, m, params, x = seq(0, 100)) {
  sim <- expand.grid(x = x, param = names(params), lvl = seq(0, 1, 0.1))
  par_val <- matrix(0.5, ncol = length(params), nrow = nrow(sim))
  colnames(par_val) <- names(params)
  sim <- cbind(sim, par_val)
  for (coln in names(params)) {
    sim[sim$param == coln, coln] <-
      sim[sim$param == coln, "lvl"] * (diff(params[[coln]]) + params[[coln]][1])
    sim[sim$param == coln, setdiff(names(params), coln)] <-
      data.frame(lapply(params[setdiff(names(params), coln)], mean))
  }
  sim$stocks <- do.call(y, sim[, c("x", names(params))])
  sim$productivity <- do.call(p, sim[, c("x", names(params))])
  sim$mortality <- do.call(m, sim[, c("x", names(params))])

  param_name <- data.frame(
    param = names(params),
    plabel = paste0(
      names(params), " in [",
      lapply(params, function(x) x[1]),
      ";",
      lapply(params, function(x) x[2]),
      "]"
    )
  )
  sim <- sim |>
    pivot_longer(cols = c("stocks", "productivity", "mortality")) |>
    merge(param_name)
  sim$name <- factor(sim$name, c("stocks", "productivity", "mortality"))
  ggplot(sim, aes(x, value, col = lvl, group = lvl)) +
    geom_line() +
    labs(x = NULL, y = NULL, col = "Parameter\nvalue") +
    scale_color_viridis_c(
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Medium", "High")
    ) +
    facet_grid(name ~ plabel, scales = "free") +
    theme_minimal()
}


## test
y <- function(x, theta, lambda, tau, kappa, omega, alpha, beta, y0) {
  y0 + (theta - y0) / (1 + kappa * exp(-lambda * (x - tau))) + alpha /
    beta * exp(-beta * x)
}
p <- function(x, theta, lambda, tau, kappa, omega, alpha, beta, y0) {
  (theta - y0) * kappa * lambda * exp(-lambda * (x - tau)) /
    (1 + kappa * exp(-lambda * (x - tau)))^2 +
    omega * y(x, theta, lambda, tau, kappa, omega, alpha, beta, y0)
}
m <- function(x, theta, lambda, tau, kappa, omega, alpha, beta, y0) {
  omega * y(x, theta, lambda, tau, kappa, omega, alpha, beta, y0) +
    alpha * exp(-beta * x)
}
params <- list(
  theta = c(100, 400),
  lambda = c(0.1, 1),
  tau = c(5, 50),
  kappa = c(0.1, 10),
  omega = c(0.005, 0.1),
  alpha = c(5, 200),
  beta = c(0.05, 0.5),
  y0 = c(0, 100)
)
visualisation(y, p, m, params)
