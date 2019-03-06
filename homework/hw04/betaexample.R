
## Simulate some data for illustration

n <- 100; alpha <- 3; beta <- 4
x <- rbeta(100, alpha, beta)
hist(x)

## Write a function for the negative log-likelihood

nll <- function(par, x, verbose = FALSE){
  alpha <- par[1]; beta <- par[2] # unpack
  ll <- sum(dbeta(x, alpha, beta, log = TRUE))
  if (verbose) print(c(par, -ll))
  return(-ll)
}

## Numerically minimize it

start <- c(alpha = 1, beta = 1) # starting values
eps <- 1e-10 # small value for lower bounds

op <- optim(par = start, fn = nll,
            lower = rep(eps, 2),
            x = x, verbose = TRUE)
op # This is a list; extract elements using $
mle <- op$par

## Examine the nll at the min

alpha.test <- cbind(seq(mle[1]/2, mle[1]*2, length = 100), mle[2])
nll.alpha <- apply(alpha.test, 1, nll, x = x, verbose = FALSE)
plot(alpha.test[,1], nll.alpha, type = "l")

beta.test <- cbind(mle[1], seq(mle[2]/2, mle[2]*2, length = 100))
nll.beta <- apply(beta.test, 1, nll, x = x, verbose = FALSE)
plot(beta.test[,2], nll.beta, type = "l")

## Do it again, this time estimating the Hessian

op <- optim(par = start, fn = nll,
            lower = rep(eps, 2), hessian = TRUE,
            x = x, verbose = TRUE)
mle <- op$par
J <- solve(op$hessian) # no negative - already working with negative ll
se.hat <- sqrt(diag(J))

lower <- mle - 2*se.hat
upper <- mle + 2*se.hat


