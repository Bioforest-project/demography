integral_stpp <- function(x,
                          delta = 0.2,
                          tau = 10,
                          alpha = 0.1,
                          lambda = 0.05,
                          omega = 0.02,
                          y0 = 10) {
  (delta * exp(2) * tau) / 4 -
    delta * exp(2 * (1 - x / tau)) * (x^2 / (2 * tau) + x / 2 + tau / 4)
}

integral_stpm <- function(x,
                          delta = 0.2,
                          tau = 10,
                          alpha = 0.1,
                          lambda = 0.05,
                          omega = 0.02,
                          y0 = 10) {
  alpha / lambda * (1 - exp(-lambda * x))
}

stocks <- function(x,
                   delta = 0.2,
                   tau = 10,
                   alpha = 0.1,
                   lambda = 0.05,
                   omega = 0.02,
                   y0 = 10) {
  y0 + integral_stpp(x, delta, tau, alpha, lambda, omega, y0) -
    integral_stpm(x, delta, tau, alpha, lambda, omega, y0)
}

prod <- function(x,
                 delta = 0.2,
                 tau = 10,
                 alpha = 0.1,
                 lambda = 0.05,
                 omega = 0.02,
                 y0 = 10) {
  delta * (x / tau)^2 * exp(2 * (1 - x / tau)) +
    omega * stocks(x, delta, tau, alpha, lambda, omega, y0)
}

mort <- function(x,
                 delta = 0.2,
                 tau = 10,
                 alpha = 0.1,
                 lambda = 0.05,
                 omega = 0.02,
                 y0 = 10) {
  alpha * exp(-lambda * x) +
    omega * stocks(x, delta, tau, alpha, lambda, omega, y0)
}

integral_stocks <- function(x,
                            delta = 0.2,
                            tau = 10,
                            alpha = 0.1,
                            lambda = 0.05,
                            omega = 0.02,
                            y0 = 10) {
  delta * exp(2 * (1 - x / tau)) * tau / 4 *
    ((x^2 / (tau) + 2 * x + tau * 3 / 2)) +
    x * (y0 + delta * exp(2) * tau / 4 - alpha / lambda) -
    exp(2) / 8 * delta * tau^2 -
    alpha / lambda^2 * exp(-lambda * x)
}

delta_prod <- function(x,
                       delta = 0.2,
                       tau = 10,
                       alpha = 0.1,
                       lambda = 0.05,
                       omega = 0.02,
                       y0 = 10,
                       k = 2) {
  (integral_stpp(x + k) - integral_stpp(x) +
     omega * (integral_stocks(x + k) - integral_stocks(x))) / k
}

delta_mort <- function(x,
                       delta = 0.2,
                       tau = 10,
                       alpha = 0.1,
                       lambda = 0.05,
                       omega = 0.02,
                       y0 = 10,
                       k = 2) {
  (integral_stpm(x + k) - integral_stpm(x) +
     omega * (integral_stocks(x + k) - integral_stocks(x))) / k
}
