# Suppose p(theta | y) = N(0, 1)
# We want to estimate P(theta > 3 | y)

# Monte Carlo integration will be a poor approximation to this
# because only a small portion of the samples from the 
# N(0, 1) will be above 3 (and the relative error about which
# side of 3 they fall will be important).

# Let g ~ N(mu, 1)

# Thus, we have:
# q(theta | y) = exp(-theta^2/2)
# g(theta) = (2 * pi)^(-1/2) * exp(-(theta - mu)^2/2)
# w(theta) = q(theta |y) / g(theta) = sqrt(2 * pi) * exp(mu*(mu - 2*theta)/2)

# We'll construct confidence bands for the mean and 
# identify the value of mu with the smallest standard error.

set.seed(1)
B = 100000
# take sample from N(0, 1) (will shift by mu to get a sample from g)
x = rnorm(B)

# x is a sample from a N(0, 1) distribution
# loc is the mean of g(theta), i.e., g ~ N(loc, 1) distribution
imp <- function(x, loc = 0) {
  B = length(x)
  theta = x + loc # create samples from g by shifting N(0, 1) sample by loc
  # compute importance ratio
  w = sqrt(2 * pi) * exp(loc*(loc - 2*theta)/2)
  t = (theta > 3) * w
  
  muhat <- sum(t)/sum(w) # estimate P(theta > 3 | y)
  # estimate variance of estimated mean
  v = 1/B * (var(t) + muhat^2 * var(w) - 2 * muhat * cov(t, w))
  # compute lower and upper bound of confidence interval
  lb = muhat - qnorm(0.975) * sqrt(v)
  ub = muhat + qnorm(0.975) * sqrt(v)
  # compute width of confidence interval
  width = ub - lb
  list(loc = loc, v = v, muhat = muhat, lb = lb, ub = ub, width = width, 
       theta = theta, w = w, t = t)
}

# do importance sampling for different values of mu
m = seq(0, 4, len = 101)
out = vector("list", length(m))
for (i in 1:length(m)) {
  out[[i]] = imp(x, m[i])
}

# plot estimates and confidence bands of P(theta > 3 | y) 
# for different values of mu

# extract and vectorize elements from out list
muhat = sapply(out, getElement, "muhat")
lb = sapply(out, getElement, "lb")
ub = sapply(out, getElement, "ub")
width = sapply(out, getElement, "width")
wmin = which.min(width) # element list with narrowest CI width

plot(m, muhat, type = "l", ylim = range(c(lb, ub)),
     ylab = expression(P(theta > 3*" | "*y)), xlab = expression(mu))
lines(m, lb, lty = 2)
lines(m, ub, lty = 2)
abline(v = m[wmin], col = "blue") # location of narrowest CI
abline(h = (1 - pnorm(3)), col = "orange") # true probability
legend("topright", 
       legend = c("estimated probability", "95% confidence bands", 
                  "narrowest CI", "true probability"), 
       lty = c(1, 2, 1, 1), col = c("black", "black", "blue", "orange"))
title("Importance sampling estimates\nas a function of proposal distribution")

# histogram of log weights
# the histogram is symmetric, which means the sampler should be okay
hist(log(out[[wmin]]$w), xlab = "log(w)", main = "")

# plot of t(theta) versus theta
# sequence of values over range of sampled theta values
# for sampler with narrowest CI
theta = out[[wmin]]$theta
x = seq(min(theta), max(theta), len = 1000)
loc = out[[wmin]]$loc
# This plot has a uniformly bounded spike, so the sampler
# should work reasonably well
plot(x,  (x > 3) * sqrt(2 * pi) * exp(loc*(loc - 2*x)/2), xlab = "theta",
     ylab = "t(theta)", type = "l")

