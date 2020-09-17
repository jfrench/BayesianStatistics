# An example of finding the posterior mode
# (maximum a posteriori estimation)

# Data distribution: y_1, ..., y_n | mu, sigma^2
# ~ i.i.d. N(mu, sigma^2)
# Prior distribution: mu ~ U(10, 15)
# sigma^2 ~ N(0.5, 0.25) * I(sigma^2 > 0)
# Essentially a truncated normal with a
# different scaling constant

# Generate some synthetic data
# Normal with mean 11 and sd = 0.47
set.seed(7)
y = rnorm(10, mean = 11, sd = 0.47)

# create function for negative log unnormalized posterior
# takes the vector theta = (mu, sigma^2)
# and y = data.
# Note that because of the priors, mu must
# be between 10 and 15 and sigma^2 > 0.
# This isn't explicitly captured in the function
# We multiply the log posterior by -1 so that instead of
# trying to maximize the log unnormalized posterior, we minimize
# the negative log unnormalized posterior
nlup = function(theta, y) {
  mu = theta[1]
  sigma = sqrt(theta[2])
  obj = sum(dnorm(y, mean = mu, sd = sigma, log = TRUE)) +
    dunif(mu, 10, 15, log = TRUE) +
    dnorm(sigma^2, mean = 0.5, sd = 0.5, log = TRUE)
  -obj
}

# optim performs multi-dimensional optimization
# par - vector of starting values
# f - function to MINIMIZE. The first argument
# must be the argument you want to optimize over
# lower - the constraints on the lower bound
# upper - the constraints on the upper bound
# method - the optimization method. "L-BFGS-B"
# allows you to specify constraints
# control - a list of optional tuning parameters.
# The remaining arguments are the arguments that
# must be supplied to f.
optim(par = c(12.5, 1),
      f = nlup,
      lower = c(10.0001, 0.05),
      upper = c(14.9999, 2),
      method = "L-BFGS-B",
      y = y)

# nlminb is an alternative optimizer
# the purpose of each argument should be
# straightforward
(map = nlminb(start = c(12.5, 1),
              objective = nlup,
              lower = c(10.0001, 0.0001),
              upper = c(14.9999, 10),
              y = y))

# create sequence of values for mu and sigmasq
mymu = seq(10.8, 11.2, length = 200)
mysigmasq = seq(0.15, 0.45, length = 200)
# create grid
mytheta = expand.grid(mymu, mysigmasq)
# for each row of mytheta, plug it into nlup function
z = apply(mytheta, 1, function(theta) {
  nlup(theta, y = y)
})
# convert z to matrix for plotting
zmat = matrix(z, nrow = length(mymu))

# create heat map of objective surface
image(mymu, mysigmasq, zmat,
      col = hcl.colors(64, "YlOrRd", rev = TRUE),
      xlab = expression(mu), ylab = expression(sigma^2))
# add contours
contour(mymu, mysigmasq, zmat, add = TRUE)
# place point for posterior mode
points(map$par[1], map$par[2], pch = 20)
title("log posterior density (unnormalized)")
