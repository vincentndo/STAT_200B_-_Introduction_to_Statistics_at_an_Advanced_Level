
## Generate data for example

r <- function(x){ # A test function -- try different ones
  (10*exp(-x^2) + 5*cos(2*pi*x^2))
}

n <- 100
x <- runif(n)
y <- r(x) + rnorm(n, sd = 2)

par(mar = c(5, 4, 1, 1))
plot(x, y, pch = 16, cex = 0.5)
xseq <- seq(0, 1, length = 100)
lines(xseq, r(xseq), lwd = 2)
legend(locator(1), lty = 1, bty = "n", lwd = 2, legend = "True function")
dev.print(pdf, file = "example1.pdf", height = 5, width = 7)

## k-nearest neighbors

knn <- function(x, y, xseq, k){
  require(fields)
  dmat <- rdist(x, xseq)
  indices <- apply(dmat, 2, order)[1:k,]
  return(apply(indices, 2, function(i){mean(y[i])}))
}

plot(x, y, pch = 16, cex = 0.5)
lines(xseq, r(xseq), lwd = 2)
lines(xseq, knn(x, y, xseq, k = 20), col = 2, lwd = 2)
lines(xseq, knn(x, y, xseq, k = 5), col = 4, lty = 2, lwd = 2)
legend(locator(1), lwd = 2, lty = c(1, 1, 2), col = c(1, 2, 4), bty = "n",
       legend = c("True function", "k = 20", "k = 5"), y.intersp = 2)
dev.print(pdf, file = "example2.pdf", height = 5, width = 7)

## N-W kernel estimator, normal kernel

nwkernel <- function(x, y, xseq, h){
  require(fields)
  dmat <- rdist(x, xseq)
  K <- dnorm(dmat/h)
  rhat <- sapply(1:length(xseq), function(j){
    sum(K[,j]/sum(K[,j]) * y)
  })
  return(rhat)
}

plot(x, y, pch = 16, cex = 0.5)
lines(xseq, r(xseq), lwd = 2)
lines(xseq, nwkernel(x, y, xseq, h = 0.1), lwd = 2, col = 2)
lines(xseq, nwkernel(x, y, xseq, h = 0.01), lwd = 2, col = 4, lty = 2)
legend(locator(1), lwd = 2, lty = c(1, 1, 2), col = c(1, 2, 4), bty = "n",
       legend = c("True function", "h = 0.1", "h = 0.01"), y.intersp = 2)
dev.print(pdf, file = "example3.pdf", height = 5, width = 7)

nwk.risk <- function(h, x, y){
  require(fields)
  dmat <- rdist(x)
  K <- dnorm(dmat/h)
  rhat <- sapply(1:length(x), function(j){
    sum(K[,j]/sum(K[,j]) * y)
  })
  sum((y-rhat)^2 / (1 - dnorm(0)/apply(K, 1, sum))^2)
}

hseq <- seq(0.001, 0.1, length = 100)
risk <- sapply(hseq, nwk.risk, x = x, y = y)
plot(hseq, risk, type = "l")

hopt <- optimize(nwk.risk, lower = 0.001, upper = 0.5, x = x, y = y)$min

plot(x, y, pch = 16, cex = 0.5)
lines(xseq, r(xseq), lwd = 2)
lines(xseq, nwkernel(x, y, xseq, h = hopt), lwd = 2, col = 2)
legend(locator(1), lwd = 2, lty = 1, col = 1:2, bty = "n",
       legend = c("True function", "Optimal h"), y.intersp = 2)
dev.print(pdf, file = "example4.pdf", height = 5, width = 7)

## Illustrate cosine basis

cosfuncs <- sqrt(2) * cos(outer(xseq, pi*(1:4), FUN = "*"))
cosfuncs <- cbind(rep(1, length(xseq)), cosfuncs)
matplot(xseq, cosfuncs, xlab = "x", ylab = "", type = "l", lwd = 2)
dev.print(pdf, file = "cosbasis.pdf", height = 6, width = 8)

## Illustrate Legendre polynomials

xseq <- seq(-1, 1, length = 100)
J <- 5
legpolys <- matrix(NA, nrow = length(xseq), ncol = J)
legpolys[,1] <- 1
legpolys[,2] <- xseq
for(j in 2:(J-1)) legpolys[,j+1] <- ((2*j+1)*xseq*legpolys[,j] - j*legpolys[,j-1])/(j+1)
for(j in 1:J) legpolys[,j] <- sqrt((2*j+1)/2) * legpolys[,j]
matplot(xseq, legpolys, type = "l", lwd = 2, xlab = "x", ylab = "")
dev.print(pdf, file = "legpolys.pdf", height = 6, width = 8)

## Estimate the test function using orthogonal series regression - use cosine basis

osreg <- function(x, y, xseq, J){
  Phi <- sqrt(2) * cos(outer(x, pi*(1:(J-1)), FUN = "*"))
  Phi <- cbind(rep(1, length(xseq)), Phi)
  beta.hat <- lm(y~Phi-1)$coef
  Phi.pred <- sqrt(2) * cos(outer(xseq, pi*(1:(J-1)), FUN = "*"))
  Phi.pred <- cbind(rep(1, length(xseq)), Phi.pred) 
  return(Phi.pred %*% beta.hat)
}

plot(x, y, pch = 16, cex = 0.5)
lines(xseq, r(xseq), lwd = 2)
lines(xseq, osreg(x, y, xseq, J = 4), col = 2, lwd = 2)
lines(xseq, osreg(x, y, xseq, J = 10), col = 4, lwd = 2, lty = 2)
legend(locator(1), lwd = 2, lty = c(1, 1, 2), col = c(1, 2, 4), bty = "n",
       legend = c("True function", "J = 4", "J = 10"), y.intersp = 2)
dev.print(pdf, file = "example5.pdf", height = 5, width = 7)

osreg.risk <- function(J, x, y){
  Phi <- sqrt(2) * cos(outer(x, pi*(1:(J-1)), FUN = "*"))
  Phi <- cbind(rep(1, length(xseq)), Phi)
  mod <- lm(y~Phi-1)
  u <- lm.influence(mod)$hat
  return(sum((mod$res/(1-u))^2))
}

Jvals <- 1:15
rhat <- sapply(Jvals, osreg.risk, x = x, y = y)
plot(Jvals, rhat, type = "h")
Jopt <- Jvals[rhat == min(rhat)]
plot(x, y, pch = 16, cex = 0.5)
lines(xseq, r(xseq), lwd = 2)
lines(xseq, osreg(x, y, xseq, J = Jopt), col = 2, lwd = 2)
legend(locator(1), lwd = 2, lty = 1, col = 1:2, bty = "n",
       legend = c("True function", "Optimal J"), y.intersp = 2)
dev.print(pdf, file = "example6.pdf", height = 5, width = 7)
