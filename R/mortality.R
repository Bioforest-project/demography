mortality <- function(x,
                      tau = 10,
                      lambda = 0.05,
                      delta = 1,
                      thetainf = 0.5,
                      theta0 = 0.4,
                      omega = 0.02,
                      ba0 = 10) {
  k <- ba0 - thetainf / omega - ((theta0 - thetainf) * lambda) /
    (omega * lambda - 1) -
    (16 * exp(2) * delta * (thetainf - theta0)) / (tau^2 * (2 * omega - tau)^3)
  stocks <- thetainf / omega + k * exp(-omega * x) +
    ((theta0 - thetainf) * lambda) / (omega * lambda - 1) * exp(-lambda * x) +
    (exp(1 - x / tau))^2 * (2 * delta * (thetainf - theta0)) /
      (tau^2 * (2 * omega - tau)) * (x^2 - x * 4 / (2 * omega - tau) +
                                       8 / (2 * omega - tau)^2)
  return(omega * stocks)
}
