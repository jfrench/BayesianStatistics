# Suppose q(theta | y) is proportional to exp(-|theta - 5|).
# This is the kernel of a double exponential distribution.

# We want to determine the second moment of the 
# target distribution, i.e., E(h(theta) | y) 
# with h(theta) = theta^2.

set.seed(2)

# unnormalized target density
qtarget = function(theta) exp(-abs(theta - 5))

# Estimate second moment analytically 
# (1/2 is the normalizing constant)
f = function(theta) qtarget(theta)/2*theta^2
(truth = integrate(f, lower = -5, upper = 15)$value)

# plot qtarget
x = seq(-5, 15, len = 1001)
plot(x, qtarget(x), type = "l", 
     xlab = expression(x), 
     ylab = expression(q(theta*"|"*y)))

# For our proposal distribution (g(theta)), we'll use a N(5, 2^2)
# compare q(theta | y) to g(theta)
lines(x, dnorm(x, mean = 5, sd = 2), col = "orange")

# h(theta)q(theta|y)/g(theta) is somewhat flat over the
# most likely part of the domain, but has unbounded spikes in 
# the tails, which is a bad sign
plot(x, x^2*qtarget(x)/dnorm(x, mean = 5, sd = 2), 
     col = "blue", type = "l", ylim = c(0, 2000), 
     xlab = expression(theta), ylab = "hq/g")

# do importance sampling using g
B = 1e6 # Determine number of samples
theta = rnorm(B, mean = 5, sd = 2) # draw theta sample from g
w = qtarget(theta)/dnorm(theta, mean = 5, sd = 2) # calculate importance weights
t = theta^2 * w
(muhat <- sum(t)/sum(w)) # estimate second moment
muhat1 <- muhat # store for comparison later
# estimated variance of estimate
(v = 1/B*(var(t) + muhat^2*var(w) - 2 * muhat * cov(t, w)))
v1 <- v # store for comparison later

# histogram of log weights
# this is positively skewed, which is a bad sign
hist(log(w), xlab = "log(w)")

# Another try with a different proposal density 
# (a scaled t w/ 2 df, and mean 5)
source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")

# compare q with new g
plot(x, qtarget(x), type = "l", xlab = expression(theta), 
     ylab = expression(q(theta*"|"*y)))
lines(x, dst(x, mean = 5, df = 2), col = "orange")

# h(theta)q(theta|y)/g(theta) plot should be flat
# compare this for new g and old g
# the plot has a uniform bound for the g, but not the old g
# this suggests our new sampler should work better than our previous
# sampler
plot(x, x^2*qtarget(x)/dst(x, mean = 5, df = 2), 
     col = "blue", type = "l", ylim = c(0, 200), 
     xlab = expression(theta),
     ylab = "hq/g")
lines(x, x^2*qtarget(x)/dnorm(x, mean = 5, sd = 2), col = "orange") 
legend("topright", legend = c("N(5, 2^2)", expression(t[2](5))), 
       col = c("orange", "blue"), lty = 1)

# Determine number of samples
theta = rst(B, df = 2, mean = 5) # simulate from g
w = qtarget(theta)/dst(theta, mean = 5, df = 2) # calculate importance weights
t = theta^2 * w
(muhat <- sum(t)/sum(w)) # estimate second moment
muhat2 <- muhat # store for later comparison
# estimated variance of estimate
(v = 1/B*(var(t) + muhat^2*var(w) - 2 * muhat * cov(t, w)))
v2 <- v # store for later comparison

# histogram of log weights
# this histogram is negatively skewed, which suggests
# our sampler should work reasonably well
hist(log(w), xlab = "log(w)")

# compare truth with estimate and estimated standard error
c(truth, muhat1, muhat2)
c(v1, v2)
c(truth - muhat1, truth - muhat2)

# compare target with proposal densities on common scale
plot(x, qtarget(x), type = "l", xlab = expression(theta), 
     ylab = expression(q(theta*"|"*y)))
dn = dnorm(x, mean = 5, sd = 2)
lines(x, dn/max(dn), col = "orange")
dt2 = dst(x, mean = 5, df = 2)
lines(x, dt2/max(dt2), col = "blue")
legend("topleft", legend = c("target", "normal", "t"),
	lty = 1, col = c("black", "orange", "blue"))

