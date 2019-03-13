
n <- 100 # Sample size in each simulation
B <- 100000 # Number of simulations

nullprop <- 0.99 # Proportion of times the null is actually true (mu = 0)
nullcount <- floor(nullprop * B)
mu <- c(rep(0, nullcount), runif(B-nullcount, -1, 1)) # Can replace the uniforms with something else

pvals <- sapply(1:B, function(i){
  x <- rnorm(n, mean = mu[i], sd = 1)
  w <- mean(x)/sqrt(var(x)/n)
  return(2*pnorm(-abs(w)))
})

check <- which(pvals > 0.01 & pvals < 0.05) # Indices for tests to look at
print(mean(check <= nullcount)) # Condition satisfied means H0 was true
