### Two examples of numerical integration

### estimating the mean of a distribution for a N(0, 1)
# stochastically
mean(rnorm(10))
mean(rnorm(100))
mean(rnorm(1000))
mean(rnorm(10000))
mean(rnorm(100000))

# deterministically
e = function(x) x * dnorm(x)
integrate(e, lower = -10, upper = 10)