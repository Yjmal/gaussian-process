N <- 1000

xSim <- matrix(runif(3*N, -pi, pi), N, 3)

ishigamiFun <- function(x){
  sin(x[, 1]) + 7 * sin(x[, 2])^2 + 0.1 * x[, 3]^4 * sin(x[, 1])
}

ySim <- ishigamiFun(xSim)

par(mfrow = c(1, 3))

for (i in 1:3){
  xi <- xSim[, i]
  plot(xi, ySim)
  ss <- loess(y ~ xi, data = data.frame(y = ySim, xi = xi)) # estimation of the conditional expectation
  t <- seq(from = min(xi), to = max(xi), length = 200)
  lines(t, predict(ss, t), col = "red", lwd = 2)
}

