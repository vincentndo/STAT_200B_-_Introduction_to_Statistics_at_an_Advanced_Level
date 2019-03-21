
## Specify the prior distribution

m <- .40  # Prior mean
v <- .00001 # Prior variance

beta <- (m^2*(1-m) - m*v)/v
alpha <- beta*m/(1-m)

## Plot the prior density

x <- seq(0.01, 0.99, length = 99)
prior <- dbeta(x, shape1 = alpha, shape2 = beta) # the prior pdf
plot(x, prior, type = "l", xlab = expression(theta),
     ylab = expression(p(theta)), main = "Prior density")

## The data

n <- 3547 # Number of people who took the survey
y <- round(0.49*n, 0) # Gallup reported 49% approval

## Update the parameters to get the posterior distribution

alpha.star <- y + alpha
beta.star <- n - y + beta
posterior <- dbeta(x, shape1 = alpha.star,
                   shape2 = beta.star)

## Plot the prior and posterior densities

plot(x, prior, type = "l", xlab = expression(theta), ylab = "",
     ylim = range(c(prior, posterior)))
lines(x, posterior, col = 2)
legend("topleft", lty = rep(1, 2), col = 1:2,
       legend = c("Prior", "Posterior"), bty = "n")
points(m, 0)
points(y/n, 0, col = 3)
points(alpha.star/(alpha.star+beta.star), 0, col = 2)

## Exercise: calculate the equal tail and HPD credible intervals
## (They will be close, since the posterior is nearly symmetric)
