### Sampling from a Beta(3, 3) assuming we know the
### unnormalized density.

### p(theta | y) propto theta^2 * (1 - theta)^2 = q(theta)

### Since this has finite support, our bounding distribution
### can be Uniform(0,1)
### with the appropriate bounding constant M.

# Taking the derivative of q(theta) with respect to theta,
# we can determine that the function has a maximum at theta = 1/2.
# Also, this is apparent by symmetry.  Thus, the maximum is q(1/2) = 0.0625.
# Thus, q <= g*M with M = 0.0625.

# define unnormalized density
qtarget = function(theta) {
	theta^2*(1 - theta)^2 * (theta > 0) * (theta < 1)
}

# define bounding function
gM = function(theta) {
	0.04 * dnorm(theta, mean = 0.5, sd = 1/4)
}

# plot q and bounding function
theta = seq(0, 1, len = 1000)
plot(theta, qtarget(theta), type = "l", xlab = expression(theta),
     ylab = expression(q(theta)))
lines(theta, gM(theta), col = "blue")
abline(v = 1/2)

B = 100000
mytheta = numeric(B)

i = 0 # the samples accepted
while (i < B) {
  x = rnorm(1, mean = 0.5, sd = 1/4) # sample from g distribution

  # accept x with probability q(x)/gM(x)
	if (runif(1) <= qtarget(x)/gM(x)) {
		i = i + 1
		mytheta[i] = x
	}
}

dmytheta = density(mytheta)
dtruth = dbeta(theta, 3, 3)

plot(theta, dtruth, col = "orange", type = "l")
lines(dmytheta)
legend("topleft", legend = c("approximation", "truth"), col = c("black", "orange"),
	lwd = c(1, 1))