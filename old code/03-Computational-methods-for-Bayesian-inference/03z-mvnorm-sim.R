# How to generate a sample from a multivariate normal distribution

# Assume theta is a 3x1 vector of observations.

# Let's generate an arbitrary mean vector and covariance matrix for 
# theta (where the covariance matrix must be positive definite.)

set.seed(45)
mu0 <- runif(3) # true mean
A <- matrix(runif(9), nrow = 3) # matrix of random values
Sigma <- A %*% t(A) # create pos. definite covariance matrix
Sigma

# Find square root of Sigma
Lambda.sqrt <- t(chol(Sigma))

# Verify that this has the required properties
Lambda.sqrt %*% t(Lambda.sqrt) - Sigma

nsim <- 10000

# sample values from a N(0, 1), organize in matrix
Z <- matrix(rnorm(length(mu0) * nsim), nrow = length(mu0), ncol = nsim)

# Observations from N(mu0, Sigma)
U <- mu0 + Lambda.sqrt %*% Z

# Verify that sample means and covariance matrix are close to true values
rowMeans(U) - mu0
var(t(U)) - Sigma
