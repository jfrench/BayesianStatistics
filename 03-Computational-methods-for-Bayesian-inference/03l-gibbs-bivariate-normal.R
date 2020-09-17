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

#observed data
y1 = 0
y2 = 0

#create matrix to store sims
thetapost = matrix(0, nrow = B, ncol = 2)
theta = c(-2.5, -2.5)

for (i in 1:B) {
	# determine full conditional mean
	m1 = y1 + rho*(theta[2] - y2)
	# simulate from distribution from full conditional for theta1
	theta[1] = rnorm(1, m1, sigma)

	# do the same thing for theta2
	m2 = y2 + rho*(theta[1] - y1)
	theta[2] = rnorm(1, m2, sigma)

	# save sample
	thetapost[i, ] = theta
}

#plot observations for posterior
plot(thetapost, pch = ".", xlab = expression(theta[1]),
     ylab = expression(theta[2]))
title("Samples from Gibbs sampler")

# Create three more chains
thetapost2 = thetapost3 = thetapost4 =  matrix(0, nrow = B, ncol = 2)
theta_b = c(-2.5, 2.5)
theta_c = c(2.5, -2.5)
theta_d = c(2.5, 2.5)

for (i in 1:B) {
	#determine full conditional mean,
	#simulate from distribution for theta1, theta2
	m1b = y1 + rho*(theta_b[2] - y2)
	m1c = y1 + rho*(theta_c[2] - y2)
	m1d = y1 + rho*(theta_d[2] - y2)
	theta_b[1] = rnorm(1, m1b, sigma)
	theta_c[1] = rnorm(1, m1c, sigma)
	theta_d[1] = rnorm(1, m1d, sigma)

	m2b = y2 + rho*(theta_b[1] - y1)
	m2c = y2 + rho*(theta_c[1] - y1)
	m2d = y2 + rho*(theta_d[1] - y1)
	theta_b[2] = rnorm(1, m2b, sigma)
	theta_c[2] = rnorm(1, m2c, sigma)
	theta_d[2] = rnorm(1, m2d, sigma)

	#store simulations
	thetapost2[i, ] = theta_b
	thetapost3[i, ] = theta_c
	thetapost4[i, ] = theta_d
}

#plot samples from each chain
plot(thetapost, pch = ".",
     xlab = expression(theta[1]), ylab = expression(theta[2]))
points(thetapost2, pch = ".", col = "orange")
points(thetapost3, pch = ".", col = "blue")
points(thetapost4, pch = ".", col = "grey")
legend("topleft", col = c("black", "orange", "blue", "grey"),
	pch = 20, legend = c("Chain 1", "Chain 2", "Chain 3", "Chain 4"))
title("Samples from Gibbs sampler")

# initial values of chains
theta = c(-2.5, -2.5)
theta2 = c(-2.5, 2.5)
theta3 = c(2.5, -2.5)
theta4 = c(2.5, 2.5)

#plot path of first 10 cycles of all chains
#plot.mcmc.path takes
# x (the mcmc chain)
# x0 the starting values for the chain, if desired
# other arguments to pass to the plot
# if add = TRUE, then the lines are added to the existing plot
plot.mcmc.path(thetapost[1:10,], x0 = theta,
     xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
     xlab = expression(theta[1]), ylab = expression(theta[2]))
plot.mcmc.path(thetapost2[1:10,], x0 = theta2, col = "red",
              add = TRUE)
plot.mcmc.path(thetapost3[1:10,], x0 = theta3, col = "blue",
               add = TRUE)
plot.mcmc.path(thetapost4[1:10,], x0 = theta4, col = "green",
               add = TRUE)

# plot cyles of all chains
# combine all chains
allpost = rbind(thetapost, thetapost2, thetapost3, thetapost4)
#plot first all cycles of all chains
plot(rbind(theta, thetapost), type = "l",
     xlim = range(allpost[,1]), ylim = range(allpost[,2]),
     xlab = expression(theta[1]), ylab = expression(theta[2]))
lines(rbind(theta2, thetapost2), col = "red")
lines(rbind(theta3, thetapost3), col = "blue")
lines(rbind(theta4, thetapost4), col = "green")
title("Cycles of 4 different chains")
