### Example of Gibbs sampler using bivariate normal
### Taken Section 11.3 of Bayesian Data Analysis, 3e by Gelman et al. (2013)

# load useful distributions
source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")

# Data distribution: (y1, y2) | theta = (theta1, theta2) ~
# N((theta1, theta2), ((1, rho), (rho, 1))) (rho known)
#
# Prior distribution:
# p(theta1) = p(theta2) = 1.
#
# Posterior distribution:
# p(theta1, theta2 | y) = N((y1, y2), ((1, rho), (rho, 1))).

# Full conditional distributions:
# theta1 | theta2, y = N(y1 + rho(theta2 - y2), 1 - rho^2)
# theta2 | theta1, y = N(y2 + rho(theta1 - y1), 1 - rho^2)

# set parameters
B = 1000
rho = .8
sigma = sqrt(1 - rho^2)
# observed data
y1 = 0
y2 = 0

# function to run gibbs sampler
# note that B, sigma, y1, and y2 are
# treated as global variables.
gibbs = function(theta) {
  #create matrix to store samples
  theta_sims = matrix(0, nrow = B, ncol = 2)
  # run gibbs sampler for B cycles
  for (i in 1:B) {
    # determine full conditional mean
    m1 = y1 + rho*(theta[2] - y2)
    # simulate from distribution from full conditional for theta1
    theta[1] = rnorm(1, m1, sigma)
    # do the same thing for theta2
    m2 = y2 + rho*(theta[1] - y1)
    theta[2] = rnorm(1, m2, sigma)
    # save sample
    theta_sims[i, ] = theta
  }
  return(theta_sims)
}

# run chain with an initial starting value
chain1 = gibbs(c(-2.5, 2.5))
#plot observations for posterior
plot(chain1, pch = ".", xlab = expression(theta[1]),
     ylab = expression(theta[2]))
title("Samples from Gibbs sampler")

# Create three more chains with different starting values
chain2 = gibbs(c(-2.5, 2.5))
chain3 = gibbs(c(2.5, -2.5))
chain4 = gibbs(c(2.5, 2.5))
#plot samples from each chain
plot(chain1, pch = ".",
     xlab = expression(theta[1]), ylab = expression(theta[2]))
points(chain2, pch = ".", col = "orange")
points(chain3, pch = ".", col = "blue")
points(chain4, pch = ".", col = "grey")
legend("topleft", col = c("black", "orange", "blue", "grey"),
	pch = 20, legend = c("Chain 1", "Chain 2", "Chain 3", "Chain 4"))
title("Samples from Gibbs sampler")

# initial values of chains
theta1 = c(-2.5, -2.5)
theta2 = c(-2.5, 2.5)
theta3 = c(2.5, -2.5)
theta4 = c(2.5, 2.5)

# plot path of first 10 cycles of all chains
# plot.mcmc.path takes
# x (the mcmc chain)
# x0 the starting values for the chain, if desired
# other arguments to pass to the plot
# if add = TRUE, then the lines are added to the existing plot
plot.mcmc.path(chain1[1:10,], x0 = theta1,
     xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
     xlab = expression(theta[1]),
     ylab = expression(theta[2]),
     main = "First 10 cycles of each chain")
plot.mcmc.path(chain2[1:10,], x0 = theta2, col = "red",
              add = TRUE)
plot.mcmc.path(chain3[1:10,], x0 = theta3, col = "blue",
               add = TRUE)
plot.mcmc.path(chain4[1:10,], x0 = theta4, col = "green",
               add = TRUE)

# plot cycles of all chains
# combine all chains
allpost = rbind(chain1, chain2, chain3, chain4)
#plot first all cycles of all chains
plot(rbind(theta1, chain1), type = "l",
     xlim = range(allpost[,1]), ylim = range(allpost[,2]),
     xlab = expression(theta[1]), ylab = expression(theta[2]))
lines(rbind(theta2, chain2), col = "red")
lines(rbind(theta3, chain3), col = "blue")
lines(rbind(theta4, chain4), col = "green")
title("Cycles of 4 different chains")
