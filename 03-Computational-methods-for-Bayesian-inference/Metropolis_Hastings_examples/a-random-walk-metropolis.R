library(mvtnorm) # need to sample from multivariate normal and evaluate
                 # density of multivariate normal

# Data distribution: p(y | theta) = bivariate N(y | theta, I),
# where I is the 2x2 identity matrix and y = c(0, 0) is a 2-dimensional
# vector of 0s.
# Prior distribution: the improper prior p(theta) = 1
#
# We want to sample from the joint posterior for theta.
#
# Proposal distribution:
# The proposal distribution will also be bivariate normal,
# centered at the current iteration value of theta with a
# scaled version of I for the covariance.
# i.e., theta.star | theta^(t-1) ~ N(theta^(t-1) , const * I).
#
# Since the proposal distribution is not only symmetric,
# but depends only on the distance |theta.star - theta.current|,
# we can implement the  random walk metropolis algorithm
# because the ratio r simplifies to the ratio of posterior
# densities.

#Define function to perform metropolis algorithm
mh = function(B, start, Imat, const) {
  # store iterations
  theta = matrix(0, nrow = B + 1, ncol = length(start))
  theta[1, ] = start # first row of matrix is starting values
	for (i in 2:(B + 1)) {
	  # use c to convert from matrix to vector
		theta.star = c(rmvnorm(1, theta[i - 1, ], const * Imat))
    # observed values of c(0, 0)
		r = dmvnorm(theta.star, mean = c(0, 0), sigma = Imat) /
         dmvnorm(theta[i - 1,], mean = c(0, 0), sigma = Imat)
		if (runif(1) <= min(r, 1)) {
			theta[i,] = theta.star
		} else {
			theta[i,] = theta[i - 1,]
		}
	}
	return(theta)
}

B = 1000
Imat = matrix(c(1, 0, 0, 1), nrow = 2)
mh1 = mh(B, c(0, 0), Imat, const = .2^2)
mh2 = mh(B, c(-3, 3), Imat, const = .2^2)
mh3 = mh(B, c(-3, -3), Imat, const = .2^2)
mh4 = mh(B, c(3, -3), Imat, const = .2^2)
mh5 = mh(B, c(3, 3), Imat, const = .2^2)

# Plot first 50 steps of M-H chain.  Note starting point.
plot(c(-4, 4), c(-4, 4), type = "n", xlab = "theta1", ylab = "theta2")
lines(mh1[1:51,])
lines(mh2[1:51,])
lines(mh3[1:51,])
lines(mh4[1:51,])
lines(mh5[1:51,])
points(c(-3, -3, 0, 3, 3), c(-3, 3, 0, -3, 3), pch = 20)
title("First 50 iterations of M-H chain")

# Plot all steps of M-H chain.  Note starting point.  Converges
# eventually
plot(c(-4, 4), c(-4, 4), type = "n", xlab = "theta1", ylab = "theta2")
lines(mh1)
lines(mh2)
lines(mh3)
lines(mh4)
lines(mh5)
points(c(-3, -3, 0, 3, 3), c(-3, 3, 0, -3, 3), pch = 20)
title("First 1000 iterations of M-H chain")

#plot second half of posterior sample
plot(mh1[(B/2 + 1):(B + 1), ], pch = ".", xlim = c(-4, 4), ylim = c(-4, 4),
     xlab = expression(theta[1]), ylab = expression(theta[2]))
points(mh2[(B/2 + 1):(B + 1), ], pch = ".")
points(mh3[(B/2 + 1):(B + 1), ], pch = ".")
points(mh4[(B/2 + 1):(B + 1), ], pch = ".")
points(mh5[(B/2 + 1):(B + 1), ], pch = ".")

# Try same things with a different proposal that allows for larger jumps
mh1 = mh(B, c(0, 0), Imat, const = 1)
mh2 = mh(B, c(-3, 3), Imat, const = 1)
mh3 = mh(B, c(-3, -3), Imat, const = 1)
mh4 = mh(B, c(3, -3), Imat, const = 1)
mh5 = mh(B, c(3, 3), Imat, const = 1)

# Plot first 50 steps of M-H chain.  Note starting point.
plot(c(-4, 4), c(-4, 4), type = "n",
     xlab = expression(theta[1]), ylab = expression(theta[2]))
lines(mh1[1:51,])
lines(mh2[1:51,])
lines(mh3[1:51,])
lines(mh4[1:51,])
lines(mh5[1:51,])
points(c(-3, -3, 0, 3, 3), c(-3, 3, 0, -3, 3), pch = 20)
title("First 50 iterations of M-H chain")

# Plot all steps of M-H chain.  Note starting point.  Converges
# eventually
plot(c(-4, 4), c(-4, 4), type = "n",
     xlab = expression(theta[1]), ylab = expression(theta[2]))
lines(mh1)
lines(mh2)
lines(mh3)
lines(mh4)
lines(mh5)
points(c(-3, -3, 0, 3, 3), c(-3, 3, 0, -3, 3), pch = 20)
title("Iterations of M-H chain")

#plot second half of posterior sample
plot(mh1[(B/2 + 1):(B + 1), ], pch = ".", xlim = c(-4, 4), ylim = c(-4, 4),
     xlab = expression(theta[1]), ylab = expression(theta[2]))
points(mh2[(B/2 + 1):(B + 1), ], pch = ".")
points(mh3[(B/2 + 1):(B + 1), ], pch = ".")
points(mh4[(B/2 + 1):(B + 1), ], pch = ".")
points(mh5[(B/2 + 1):(B + 1), ], pch = ".")