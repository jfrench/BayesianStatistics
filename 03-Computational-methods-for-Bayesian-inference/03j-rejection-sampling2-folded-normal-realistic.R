### Sampling from a folded normal N(0, 1).
### If y ~ N(0, 1), then |y| is a folded normal

### q(theta | y) = exp(-theta^2/2)
### this looks like a standard normal pdf, BUT it only has
### non-negative support

# We can envelope this using a chi-square distribution 2 df
qtarget = function(theta) {
	exp(-theta^2/2)
}

gM = function(theta) {
	2.4 * dchisq(theta, df = 2)
}

# plot q and gM as a function of theta
theta = seq(0, 5, len = 1000)
plot(theta, gM(theta), type = "l", col = "blue", xlab = expression(theta), ylab = "")
lines(theta, qtarget(theta))
legend("topright", legend = c("q", "gM"), lty = 1, col = c("black", "blue"))

B = 10000
mytheta = numeric(B)
i = 0
while (i < B) {
	x = rchisq(1, df = 2)
	if (runif(1) <= qtarget(x)/gM(x)) {
		i = i + 1
		mytheta[i] = x
	}
}

#The actual density is
dens = function(theta) {
	sqrt(2/pi)*exp(-theta^2/2)
}

#construct histogram of simulated values
hist(mytheta)

#plot approximate density
dmytheta = density(mytheta, from = 0, to = 5, cut = 0)
plot(dmytheta, xlab = "theta", ylab = "density", main = "")
lines(theta, dens(theta), col = "blue")
legend("topright", legend = c("approximation", "truth"),
	lwd = c(1, 1), col = c("black", "blue"))

#This doesn't look right because of boundary bias!
#Essentially, the density estimator doesn't know what to do around zero.
hist(mytheta, freq = FALSE, breaks = 100)
lines(theta, dens(theta), col = "blue")
legend("topright", legend = c("approximation", "truth"),
	lwd = c(1, 1), col = c("black", "blue"))


# Visually see that we are only accepting the values below the target density.
plot(theta, gM(theta), type = "l", col = "blue")
lines(theta, qtarget(theta))

# plot accepted/rejected points with different symbols
B = 100
mytheta = numeric(B)
i = 0
while (i < B) {
	x = rchisq(1, df = 2)
	u = runif(1, 0, gM(x))
	if (u <= qtarget(x)) 	{
		i = i + 1
		mytheta[i] = x
		points(x, u, pch = 20)
	} else 	{
		points(x, u)
	}
}
