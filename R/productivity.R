productivity <- function(x,
                         tau = 10,
                         lambda = 0.05,
                         delta = 0.2,
                         thetainf = 0.5,
                         theta0 = 0.4) {
  return(
    theta0 +
      (thetainf - theta0) * (1 - exp(-lambda * x)) +
      thetainf * delta * ((x / tau) * exp(1 - x / tau))^2
  )
}
