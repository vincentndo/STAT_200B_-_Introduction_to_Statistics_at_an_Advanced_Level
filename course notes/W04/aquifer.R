
## This example requires several additional R packages.
## They are all stored on CRAN: http://cran.r-project.org/
## To download & install, just type install.packages("nameofpackage"),
## eg. install.packages("scatterplot3d")

## Wolfcamp aquifer data (Cressie, 1984)

aquifer <- read.table("aquifer.dat", skip = 1, header = FALSE)
names(aquifer) <- c("x", "y", "z")
library(scatterplot3d) # Add-on package with 3d plotting functions
scatterplot3d(aquifer, type = "h",
              xlab = "Easting (mi)", ylab = "Northing (mi)",
              zlab = "Piezometric head (1000s ft)",
              main = "Wolfcamp Aquifer Data")
dev.print(pdf, file = "aquiferpoints.pdf", height = 8, width = 8)

## Define the negative log-likelihood function
## This is what we'll minimize

library(mvtnorm) # Add-on package with multivariate normal and t dists
nll <- function(param, Xmat, d, z, verbose = FALSE){
  beta <- param[1:3]; sigma2 <- param[4]; rho <- param[5]
  mu <- Xmat %*% beta
  Sigma <- sigma2 * exp(-d/rho)
  ll <- dmvnorm(z, mean = mu, sigma = Sigma, log = TRUE)
  if (verbose) print(c(param, -ll))
  return(-ll)
}

## Calculate the regression matrix and the matrix of distances

Xmat <- cbind(1, aquifer$x, aquifer$y)
d <- as.matrix(dist(aquifer[,1:2]))

## Choose starting values for the minimization

beta.ols <- lm(aquifer$z~Xmat-1)$coef # OLS regression; assuming iid errors
names(beta.ols) <- NULL
start <- c(beta0 = beta.ols[1], beta1 = beta.ols[2], beta2 = beta.ols[3],
           sigma2 = 1, rho = 1)

## Use optim to minimize

eps <- 1e-10 # A small value for the lower bound on sigma2, rho
op <- optim(start, fn = nll, method = "L-BFGS-B",
            lower = c(-Inf, -Inf, -Inf, eps, eps),
            hessian = TRUE, # Approximate matrix of 2nd derivs
            Xmat = Xmat, d = d, z = aquifer$z, verbose = TRUE)
op # This is a list; access elements using $

## Extract the MLEs and their estimated standard errors

mle <- op$par
J <- solve(op$hessian)
se.hat <- sqrt(diag(J))

## Examine the behavior of the nll at the min

par(mfrow = c(3, 2))
sapply(1:length(mle), function(i){
  values <- matrix(rep(mle, 100), nrow = 100, byrow = TRUE)
  values[,i] <- seq(mle[i]/2, mle[i]*2, length = 100)
  nll.seq <- apply(values, 1, nll, Xmat = Xmat, d = d, z = aquifer$z)
  plot(seq(mle[i]/2, mle[i]*2, length = 100), nll.seq, type = "l", xlab = names(mle)[i])
})

## Construct confidence intervals

lower <- mle - 2 * se.hat
upper <- mle + 2 * se.hat

## Kriging -- not part of the likelihood example,
## but I thought you might like to see it

## Define the locations at which to predict

gridn <- 50
xseq <- seq(-150, 150, length = gridn)
yseq <- seq(0, 200, length = gridn)
krige.grid <- expand.grid(x = xseq, y = yseq)

## Kriging using estimated parameters

dpred <- as.matrix(dist(rbind(aquifer[,1:2], krige.grid)))[1:85,-(1:85)]
Xpred <- cbind(1, krige.grid$x, krige.grid$y)
g <- exp(-dpred/mle["rho"])
G <- exp(-d/mle["rho"])
zpred <- Xpred %*% mle[1:3] + t(g) %*% solve(G, aquifer$z - Xmat %*% mle[1:3])

## Plotting - restrict to convex hull of the data

library(fields) # Add-on package for spatial stats
zmat <- matrix(zpred, nrow = 50)
ok <- in.poly(krige.grid, aquifer[,1:2], convex.hull = TRUE)
zmat[!ok] <- NA
drape.plot(x = xseq, y = yseq, z = zmat,
           theta = 120, phi = 30, # Viewing angles
           col = terrain.colors(60), horiz = FALSE)
