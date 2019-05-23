
## Example, page 48

B <- 10000; y <- rnorm(B)
sqrt(2*pi) * mean(sin(y)^2)

## Example, page 50

n <- 10; lambda <- 5; B <- 10000
samples <- matrix(rexp(n*B, rate = 1/lambda), nrow = B, ncol = n)
y <- apply(samples, MARGIN = 1, FUN = median)
var(y)

## Another way to do it

y <- replicate(B, median(rexp(n, rate = 1/lambda)))
var(y)

## Old Faithful data; an example of the algorithm on page 52

geyser <- read.table("http://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat", skip = 20)

## Estimate the median waiting time using plug-in estimate

x <- geyser$waiting
theta.hat <- quantile(x, probs = 0.5, type = 1)

## Use the bootstrap to get a variance estimate

# Code for generating a single T*
x.star <- sample(x, size = length(x), replace = TRUE) # Key line of code for the bootstrap!
quantile(x.star, probs = 0.5, type = 1)

# Now repeat this B times
B <- 10000
theta.boot <- replicate(B, {x.star <- sample(x, size = length(x), replace = TRUE)
                            quantile(x.star, probs = 0.5, type = 1)})

# Approximate the variance via MC integration
v.boot <- var(theta.boot)

## Confidence intervals, three ways

alpha <- 0.05

## Normal interval

z.a2 <- qnorm(alpha/2, lower.tail = FALSE) # Could approximate by 2
theta.hat + c(-1, 1) * z.a2 * sqrt(v.boot)

## Percentile interval

quantile(theta.boot, probs = c(alpha/2, 1-alpha/2), type = 1)

## Pivotal interval

quant <- quantile(theta.boot, probs = c(1-alpha/2, alpha/2),
                  type = 1, names = FALSE)
2 * theta.hat - quant


