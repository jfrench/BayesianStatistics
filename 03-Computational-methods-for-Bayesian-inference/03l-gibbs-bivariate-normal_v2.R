### Example of Gibbs sampler using bivariate normal
### Taken Section 11.3 of Bayesian Data Analysis, 3e by Gelman et al. (2013)

# load useful distributions
# devtools::install_github("jfrench/bayesutils")
library(bayesutils)

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
  #create matrix to store samples (and starting values)
  theta_sims = matrix(theta, nrow = B + 1, ncol = 2, byrow = TRUE)
  # run gibbs sampler for B cycles
  for (i in 1:B) {
    # determine full conditional mean for theta1
    m1 = y1 + rho * (theta[2] - y2)
    # simulate from full conditional distribution for theta1
    theta[1] = rnorm(1, m1, sigma)
    # determine full conditional mean for theta2
    m2 = y2 + rho * (theta[1] - y1)
    # simulate from full conditional distribution for theta1
    theta[2] = rnorm(1, m2, sigma)
    # save sample
    theta_sims[i + 1, ] = theta
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
chain2 = gibbs(c(-2.5, -2.5))
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

# create list of all chains
list_chains = list(chain1, chain2, chain3, chain4)
# plot steps of first 10 cycles of all chains
plot_mcmc_path(list_chains,
               xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
               xlab = expression(theta[1]),
               ylab = expression(theta[2]),
               main = "Steps of first 10 cycles of each chain")

# plot 100 cycles of all chains
plot_mcmc_path(list_chains,
               ncycles = 100,
               xlab = expression(theta[1]),
               ylab = expression(theta[2]),
               type = "cycle")
title("Cycles of 4 different chains")
