# library(MASS)
# library(autoimage)

### Taken from Chapter 7 of A First Course in Bayesian Statistical Methods
### by Peter Hoff.

# A sample of 22 children are given reading comprehension tests before
# and after receiving a particular instructional method.  Each student
# will then have two scores denoting the pre- and post-instructional scores

# Each row of the data set is a two-dimensional vector representing a single
# case, distributed N(mu, Sigma).

# Data distribution: bivariate normal y | mu, Sigma ~ N(mu, Sigma)

# Prior distributions:
# mu | Sigma ~ N(mu0, L0)
# Sigma ~ Inv-Wishart(nu, K0^(-1))

# The exam was designed to give average score around 50 out of 100, so
# set mu0 = c(50, 50).  Since the true mean cannot be below 0 or above
# 100, we want the prior variances to keep mu inside this range with
# high probability. We want Sigma[1,1] = Sigma[2,2] to be approximately 625, so
# that the prior probability that the mean is outside [0, 100] is only
# 0.05.  Since the two exams measure similar things, we think the
# correlation between the means is around 0.5, so
# Sigma[1,2] = Sigma[2,1] = 0.5(625) = 312.5.
# Since E(Sigma) ~ Inv-Wishart(nu0, K0^(-1)) = K0/(nu0 - d - 1),
# choosing nu0 = 4 loosely centers samples of Sigma around K0.
# So we choose K0[1,1]=K0[2,2] = 625 and the diagonal elements equal to
# 312.5.  We make the same choice for L0.

# load useful distributions
source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")

# load data
y = matrix(c(
59, 77, 43, 39, 34, 46, 32, 26, 42, 38, 38, 43, 55, 68,
67, 86, 64, 77, 45, 60, 49, 50, 72, 59, 34, 38, 70, 48,
34, 55, 50, 58, 41, 54, 52, 60, 60, 75, 34, 47, 28, 48,
35, 33), ncol = 2, byrow = TRUE)
# reformat
y = data.frame(pretest = y[,1], posttest = y[,2])

n = nrow(y)
d = ncol(y)

### Set prior parameters
mu0 = c(50, 50)
nu0 = 4
L0 = K0 = 25^2 * matrix(c(1, .5, .5, 1), nrow = 2)

### Calculate necessary summary quantities
ybar = colMeans(y)

### Determine posterior parameters
nun = nu0 + n

#set parameters
B = 10000

#create initial guess of mu, Sigma
mu = ybar
Sigma = cov(y)
# initial guess for Smu
Smu = (n - 1) * var(y)

# store matrices/arrays
mupost = matrix(0, nrow = B, ncol = d)
Sigmapost = array(0, dim = c(d, d, B))
ytildepost = matrix(0, nrow = B, ncol = d)

# Execute Gibbs sampler
for (i in 1:B) {
	#determine full conditional distribution of mean,
	#simulate from distribution
	Ln = solve(solve(L0) + n * solve(Sigma))
	mun = Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% ybar)
	mu = as.vector((rmvnorm(1, mu = mun, V = Ln)))

	#determine full conditional distribution of Sigma,
	#simulate from distribution
	Smu = (t(y) - mu) %*% t(t(y) - mu)
	Kn = K0 + Smu

	Sigma = rinvwish(1, df = nun, Sigma = Kn)[,,1]

	mupost[i, ] = mu
	Sigmapost[,,i] = Sigma

	ytildepost[i, ] = rmvnorm(1, mu = mu, V = Sigma)
}

# mean and variance of posterior samples
(mean.mupost = apply(mupost, 2, mean))
(var.mupost = apply(mupost, 2, var))
(mean.Sigmapost = apply(Sigmapost, c(1, 2), mean))
(var.Sigmapost = apply(Sigmapost, c(1, 2), var))

### quantiles of mu
apply(mupost, 2, quantile, prob = c(.01, .25, .5, .75, .99))

### plot posterior densities

#determine densities
dpretest = density(mupost[,1])
dposttest = density(mupost[,2])

plot(dpretest, type = "l", xlim = range(c(dpretest$x, dposttest$x)),
	ylim = range(c(dpretest$y, dposttest$y)), main = "")
lines(dposttest, col = "orange")
legend("topleft", legend = c("mu pretest", "mu postest"),
	col = c("black", "orange"), lwd = c(1, 1))

# Determine probability mu2 > mu1
mean(mupost[,2] > mupost[,1])

# Determine probability ytilde2 > ytilde1
mean(ytildepost[,2] > ytildepost[,1])

# compare joint density of means
d2d = MASS::kde2d(mupost[,1], mupost[,2], n = 50)
# to plot results
autoimage::pimage(d2d, col = topo.colors(64),
       xlab = "mu pretest", ylab = "mu posttest")
contour(d2d, add = TRUE)
points(mupost[1:100,], pch = 20)
