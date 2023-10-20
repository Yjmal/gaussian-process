rm(list = ls()) # to clear the environment

#### loading some packages and functions ####
# library("plot3D")
library("MASS")
library("plot3D")
source("kernFun.R")

#### Sampling from a GP  ####

## 1. Simulating sample-paths  ##
n <- 100
x <- seq(0, 1, length.out = n) # regular sequence
param <- c(1, 0.2) # covariance parameters
K1 <- expKern(x, x, param)
par(mfrow = c(1, 1))
image2D(K1, xlab = "x", ylab = "y", main = "kernel plot")
## simulating some samples using the "mvrnorm" function
samples <- mvrnorm(n = 3, mu = rep(0, n), Sigma = K1)
# ?matplot # a function to plot the samples. The samples are indexed by columns
matplot(x, t(samples), type = "l", 
        main = "sample paths", ylab = "")

## 2. Influence of a bigger lengthscale  ##
param[2] <- 3 * param[2]  
# Mulitplying by 3 the value of theta is equivalent to zooming 3 times
# on previous graphs, which can be observed visually below.
k1 <- expKern(x, x, param) 
par(mfrow = c(1,2))
image2D(k1, theta = 0, xlab = "x", ylab = "y", main = "kernel plot")
samples <- t(mvrnorm(n = 5, mu = rep(0, n), Sigma = k1))
matplot(x, samples, type = "l", main = "sample paths", ylab = "")

## 3. Matern 5/2 kernel
mat5_2Kern <- function(x, y, param){
  sigma <- param[1]
  theta <- param[2]
  d <- abs(outer(x, y, '-')) * sqrt(5) / theta
  kern <- sigma^2*(1 + d + d^2 / 3)*exp(- d)
  return(kern)
}

k2 <- mat5_2Kern(x, x, param) 
par(mfrow = c(1,1))
image2D(k2, theta = 0, xlab = "x", ylab = "y", main = "kernel plot")
samples <- t(mvrnorm(n = 5, mu = rep(0, n), Sigma = k2))
matplot(x, samples, type = "l", main = "sample paths", ylab = "")


## --------------------------- ##
## Gaussian process regression ##
## --------------------------- ##

## 5. Data construction
fun <- function(x) {
  x + sin(4*pi*x)
}
X <- seq(0.1, 1, length.out = 6)
Y <- fun(X)

## 6. Conditional mean and kernel
condMean <- function(x, X, Y, kern, param, jitter = 0) {
  K <- kern(X, X, param)
  Kinv <- solve(K) # + jitter * diag(nrow(K)))
  kxX <- kern(x, X, param)
  return(kxX %*% Kinv %*% Y)
}

condCov <- function(x, X, kern, param, jitter = 0) {
  K <- kern(X, X, param)
  Kinv <- solve(K) # + jitter *  diag(nrow(K)))
  kxX <- kern(x, X, param)
  kxx <- kern(x, x, param)
  return(kxx - kxX %*% Kinv %*% t(kxX))
}


## 7: Plot of the conditional mean / covariance
par(mfrow = c(1, 1))
x <- seq(0, 1, length.out = 250) # test points
param <- c(1, 0.1) # covariance params
kern <- mat5_2Kern   # default 
# kern <- expKern
# kern <- brownianKern   # default 
mu <- condMean(x, X, Y, kern, param, jitter = 0) # cond. mean
Sigma <- condCov(x, X, kern, param, jitter = 0) # cond. covariance

varSigma <- diag(Sigma)
varSigma <- pmax(varSigma, 0)  # to avoid numerical negative values
plot(x, fun(x), type = "l", ylim = c(-1.5, 4), ylab = "y", 
     lty = "dotted")
lines(x, mu, type = "l", col = "blue", lwd = 3)
points(X, Y, col = "black", pch = 20, cex = 3)
lines(x, mu + 1.96*sqrt(varSigma), col = "blue")
lines(x, mu - 1.96*sqrt(varSigma), col = "blue")
legend("topleft", legend = c("function f(x)", "training points",
                             "cond. mean m(x)", "95% pred. inter"),
       col = c("black", "black", "blue", "gray50"),
       lty = c(1, NaN, 1, 2), pch = c(NaN, 20, NaN, NaN))


## 8. Conditional sample-path
condSamples <- t(mvrnorm(n = 3, mu, Sigma))
matplot(x, condSamples, type = "l", ylab = "y", 
        ylim = c(-1, 4), xlim = c(0,1), col = "grey")
points(X, Y, pch = 20, cex = 2)
lines(x, fun(x), lwd = 2)
legend("topleft", legend = c("function f(x)", "training points", "cond. samples"),
       col = c("black", "black", "gray60"), lty = c(1, NaN, 2), pch = c(NaN, 20, NaN))

# remark: link with question 7
lines(x, rowMeans(condSamples), lwd = 4, col = "blue")
lines(x, mu, lwd = 4, col = "red", lty = "dotted")
lines(x, apply(condSamples, 1, quantile, 0.025), lwd = 2, col = "blue")
lines(x, apply(condSamples, 1, quantile, 0.975), lwd = 2, col = "blue")
lines(x, mu + 1.96*sqrt(varSigma), col = "red", lwd = 2, lty = "dotted")
lines(x, mu - 1.96*sqrt(varSigma), col = "red", lwd = 2, lty = "dotted")


## 9. Hint: Can you see Y in the conditional kernel?

## 10. ## We give here the code of the Brownian kernel
brownianKern <- function(x, y, param){
  outer(x, y, pmin)
}

## Bonus: GP regression with symmetry information (odd function)

kernSymm <- function(x, y, param){
  0.25 * (mat5_2Kern(x, y, param) - mat5_2Kern(-x, y, param) 
  - mat5_2Kern(x, -y, param) + mat5_2Kern(-x, -y, param) )
}
kern <- kernSymm   # default 

x <- seq(-1, 1, length.out = 250) # test points

mu <- condMean(x, X, Y, kern, param, jitter = 0) # cond. mean
Sigma <- condCov(x, X, kern, param, jitter = 0) # cond. covariance

varSigma <- diag(Sigma)
varSigma <- pmax(varSigma, 0)  # to avoid numerical negative values
plot(x, fun(x), type = "l", ylim = c(-2, 2), ylab = "y", xlim = c(-1,1))
lines(x, mu, type = "l", col = "blue")
points(X, Y, col = "black", pch = 20)
lines(x, mu + 1.96*sqrt(varSigma), col = "gray50", lty = 2)
lines(x, mu - 1.96*sqrt(varSigma), col = "gray50", lty = 2)
legend("topleft", legend = c("function f(x)", "training points",
                             "cond. mean m(x)", "95% pred. inter"),
       col = c("black", "black", "blue", "gray50"),
       lty = c(1, NaN, 1, 2), pch = c(NaN, 20, NaN, NaN))


condSamples <- t(mvrnorm(n = 5, mu, Sigma))
matplot(x, condSamples, type = "l", ylab = "y", 
        ylim = c(-2, 2), xlim = c(-1,1), lwd = 2)
points(X, Y, pch = 20, cex = 2)
lines(x, fun(x), lwd = 2)
legend("topleft", legend = c("function f(x)", "training points", "cond. samples"),
       col = c("black", "black", "gray60"), lty = c(1, NaN, 2), pch = c(NaN, 20, NaN))

