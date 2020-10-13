### Sampling from a folded normal N(0, 1).
### If y ~ N(0, 1), then |y| is a folded normal

### q(theta | y) = exp(-theta^2/2)
### this looks like a standard normal pdf, BUT it only has
### non-negative support

# We can envelope this using a standard exponential model
# with the correct multiplier

# The optimal solution is find a single intersecting
# point between the folded
# normal and our envelope (at the inflection point of the
# folded normal), which is at theta = 1.

# Setting g(theta)M = q(theta) and solving for M when
# theta = 1 results in the solution M = exp(1/2)

qtarget = function(theta) {
	exp(-theta^2/2)
}

gM = function(theta) {
	dexp(theta)*exp(1/2)
}

# plot qtarget and gM as a function of theta
theta = seq(0, 5, len = 1000)
plot(theta, qtarget(theta), ylim = c(0, qtarget(0.001)), type = "l",
     ylab = "qtarget", xlab = expression(theta))
lines(theta, gM(theta), type = "l", col = "blue")
legend("topright", legend = c("qtarget", "gM"), lty = 1, col = c("black", "blue"))

B = 10000 # number of retained samples desired
mytheta = numeric(B) # vector to store samples
i = 0 # number of retained samples
while (i < B) {
	x = rexp(1) # draw a value from proposal distribution
	# accept the value with probability based on the importance
	# ratio
	if (runif(1) <= qtarget(x)/gM(x)) {
		i = i + 1 # increment i if sample retained
		mytheta[i] = x # store sample
	}
}

#The actual density is
dens = function(theta) {
	sqrt(2/pi)*exp(-theta^2/2)
}

# plot approximate density
# This doesn't look right because of boundary bias!
# Essentially, the density estimator doesn't know what to do around zero.
# Most density estimators don't perform well when the density is non-zero at an
# end of its support, such as the exponential and half-normal densities
dmytheta = density(mytheta, from = 0, to = 5, cut = 0)
plot(dmytheta, xlab = "theta", ylab = "density", main = "")
lines(theta, dens(theta), col = "blue")
legend("topright", legend = c("approximation", "truth"),
	lwd = c(1, 1), col = c("black", "blue"))

# plot with a probability histogram instead to get a better comparison
hist(mytheta, freq = FALSE, breaks = 100)
lines(theta, dens(theta), col = "blue")
legend("topright", legend = c("approximation", "truth"),
	lwd = c(1, 1), col = c("black", "blue"))

# Visually see that we are only accepting the values below the target density.
plot(theta, gM(theta), type = "l", col = "blue", xlab = expression(theta),
     ylab = "density", ylim = c(0, 1.25))
lines(theta, qtarget(theta))
title("Accepted vs Rejected samples")

# plot accepted/rejected points with different symbols
B = 100
mytheta = numeric(B)
i = 0
while (i < B) {
  x = rexp(1)
  u = runif(1, 0, gM(x))
  if (u <= qtarget(x)) 	{
    i = i + 1
    mytheta[i] = x
    points(x, u, pch = 20)
  } else 	{
    points(x, u)
  }
}