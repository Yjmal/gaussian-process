rm(list=ls()) # to clear the environment

#### loading some packages and functions ####
library("plot3D")
library("MASS")
source("kernFun.R")

#### Example with the Exp. kernel  ####
x <- seq(0, 1, 0.01) # regular grid
param <- c(1, 0.2) # covariance parameters
k1 <- expKern(x, x, param) # computing the covariance matrix using an exp. kernel

image2D(k1, theta = 0, xlab = "x", ylab = "y") # plotting the covariance matrix
# Q: what can you observe from the covariance matrix?
?mvrnorm # using the help from RStudio

## to complete  ##
## simulating some samples using the "mvrnorm" function
# samples <- mvrnorm(...)
# ?matplot # a function to plot the samples. The samples are indexed by columns
# Q: what can you observe from the samples?
k2 <- mat5_2Kern(x,x,param)
image2D(k2, theta = 0, xlab = "x", ylab = "y")

samples <- mvrnorm(3,rep(0,101),k2)
matplot(x,t(samples))

samples2 <- mvrnorm(3,rep(0,101),k1)
matplot(x,t(samples2))
