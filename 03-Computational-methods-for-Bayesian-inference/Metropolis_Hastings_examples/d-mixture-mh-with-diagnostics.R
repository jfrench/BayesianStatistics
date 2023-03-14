library(coda)

#Example 7.2, Givens and Hoeting (2005), Estimating a mixture parameter

# Suppose we have observed 100 observations sampled
# i.i.d. from the mixture distribution: theta * N(7, .5^2) + (1 - theta) * N(10, .5^2)
# Suppose we want to estimate theta.  Note that the true theta = .7
#
# Data distribution: p(y | theta) = theta * N(7, .5^2) + (1 - theta) * N(10, .5^2)
# Prior distribution: theta ~ Beta(1, 1)
# Proposal distribution: theta_star ~ Beta(alpha, beta)

y = c(7.25325222659913, 6.85652267046824, 7.23643792894966, 7.03343611519664,
6.9186591609056, 6.65649879051228, 6.42308043084932, 7.46636287619574,
10.3497865413661, 6.93593298389149, 6.83974994639286, 10.1477534866707,
7.18844547660898, 8.79161716373787, 6.77135115622428, 9.89206349173715,
10.6292620609587, 6.17109362928208, 9.44878709751433, 7.12422462946795,
6.75066335182976, 7.42808832040163, 9.4949511197615, 6.74956775652862,
9.46445384762244, 7.27348041082583, 6.98896265672564, 7.26262394415349,
6.94244760575449, 7.28846831817204, 9.70904705672207, 10.9878054216487,
7.45111574465272, 6.97036693452533, 6.53291089305878, 6.52220443343591,
6.10163473472885, 10.2820394025033, 8.13866685075031, 9.51099560173583,
7.74154300863383, 6.14372115300404, 6.68548657458669, 10.8689484723994,
6.73064827757487, 10.5866677031468, 9.56384573435206, 6.99562383496413,
6.18576529608999, 10.5254577115642, 9.44970647261562, 9.84118730329914,
7.21312721539275, 7.83136245649827, 10.2825737379552, 6.36363038852991,
6.74285813989089, 6.98035358880607, 6.92964132433787, 6.84202550407557,
7.38016635912841, 9.78362588855716, 7.58508950152404, 9.50912675753626,
6.55132388271122, 6.88852617303433, 7.90209596243866, 7.01301915336949,
9.93871470060288, 7.51086841627442, 6.67441840075985, 6.10225392630878,
7.40059858311095, 7.30520952023732, 7.21041994910681, 10.519655486735,
10.1288125396083, 7.04575576918983, 9.8008750675778, 10.3528228749625,
7.25796001419125, 7.67456458165231, 6.52507043820612, 6.57628671699353,
6.52670898857082, 7.44576437202348, 10.5860658476161, 6.97650854570433,
6.89904378626167, 7.48642531607213, 6.39227547871176, 7.32020055238666,
9.64975498620833, 6.92112490660115, 10.4213997762887, 6.3455538999162,
6.99875637216512, 6.98618105048244, 5.96048736574766, 10.2852677946904)

# consider using two different jumping distributions:  a Beta(1,1)
# equivalent to a Uniform(0, 1), and a Beta(2, 10)
# these do not depend on the current value of theta,
# so we can implement an independence MH sampler.

#Calculate likelihood of a single value y
lik = function(y, w) {
	w * dnorm(y, 7, .5) + (1 - w) * dnorm(y, 10, .5)
}

# log of joint likelihood
log_jlik = function(y, w) {
  sum(log(lik(y, w)))
}

# plot target and proposal densities

# seq of theta values
x = seq(0, 1, len = 1001)
# plot proposal density Beta(2, 10)
plot(x, dbeta(x, 2, 10), type = "l", ylab = "density")
# Beta(1, 1) proposal density
abline(h = 1, col = "blue")
# target (scaled)
# compute log-likelihood evaluated at many values of w
loglik_x = sapply(x, function(w) sum(log(lik(y, w))))
# transform loglik_x back to original scale (with some scaling for plotting purposes)
# and plot it
lines(x, exp(loglik_x + 130), col = "orange")
# plot legend
legend("bottomright", legend = c("Beta(2, 10)", "Beta(1, 1)", "Likelihood"), col = c("black", "blue", "orange"),
       lty = 1)
title("Proposal densities vs Likelihood")
# notice that the Beta(2, 10) proposal is nowhere near the target density,
# which creates an inefficient sampler

#histogram of data
hist(y, breaks = 25, freq = FALSE)

#plot true density
theta = seq(5, 12, len = 1000)
lines(theta, lik(theta, .7))

#Derive r to understand steps
independence_chain = function(B, start, jump_parm) {
  alpha = jump_parm[1]
  beta = jump_parm[2]

  theta = numeric(B + 1)
  theta[1] = start

  for (i in 2:length(theta)) {
    theta_star = rbeta(1, alpha, beta)
    log_numr = log_jlik(y, theta_star) +
      dbeta(theta_star, 1, 1, log = TRUE) -
      dbeta(theta_star, alpha, beta, log = TRUE)
    log_denr = log_jlik(y, theta[i - 1]) +
      dbeta(theta[i - 1], 1, 1, log = TRUE) -
      dbeta(theta[i - 1], alpha, beta, log = TRUE)
    logr = log_numr - log_denr

    if (log(runif(1)) <= min(logr, 0)) {
      theta[i] = theta_star
    } else {
      theta[i] = theta[i - 1]
    }
  }
  return(theta)
}

set.seed(72)
B = 10000
ini = rbeta(1, 1, 1)
# run chains for different proposal distributions
mh1 = independence_chain(B, start = ini, jump_parm = c(1, 1))
mh2 = independence_chain(B, start = ini, jump_parm = c(2, 10))

mean(mh1)
mean(mh2)

# Histogram of values of first chain
#(after discarding the first 200 observations)
hist(mh1[201:10001], xlab = expression(theta), main = "Posterior for Beta(1,1) prior")

# Histogram of values of second chain
#(after discarding the first 200 observations)
hist(mh2[201:10001], xlab = expression(theta), main = "Posterior for Beta(2,10) prior")

# run chains from several different starting values
mh1a = independence_chain(B, start = 0, jump_parm = c(1, 1))
mh1b = independence_chain(B, start = .1, jump_parm = c(1, 1))
mh1c = independence_chain(B, start = .2, jump_parm = c(1, 1))
mh1d = independence_chain(B, start = .3, jump_parm = c(1, 1))
mh1e = independence_chain(B, start = .4, jump_parm = c(1, 1))
mh1f = independence_chain(B, start = .5, jump_parm = c(1, 1))
mh1g = independence_chain(B, start = .6, jump_parm = c(1, 1))
mh1h = independence_chain(B, start = .7, jump_parm = c(1, 1))
mh1i = independence_chain(B, start = .8, jump_parm = c(1, 1))
mh1j = independence_chain(B, start = .9, jump_parm = c(1, 1))
mh1k = independence_chain(B, start = 1, jump_parm = c(1, 1))

# create mcmc objects of each chain, discarding warmup
mc1 = mcmc(mh1a[5002:10001])
mc2 = mcmc(mh1b[5002:10001])
mc3 = mcmc(mh1c[5002:10001])
mc4 = mcmc(mh1d[5002:10001])
mc5 = mcmc(mh1e[5002:10001])
mc6 = mcmc(mh1f[5002:10001])
mc7 = mcmc(mh1g[5002:10001])
mc8 = mcmc(mh1h[5002:10001])
mc9 = mcmc(mh1i[5002:10001])
mc10 = mcmc(mh1j[5002:10001])
mc11 = mcmc(mh1k[5002:10001])

# combine chains into mcmc list
mc = mcmc.list(mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11)

# trace plot for the chains
traceplot(mc)

# acf cplot
autocorr.plot(mc)

# effective sample size
effectiveSize(mc)

# Calculate scale reduction factor
gelman.diag(mc, autoburnin = FALSE)

# Plot Gelman statistic
gelman.plot(mc, autoburnin = FALSE)

# Heidelberg & Welch
# If the halfwidth test fails, extend the chain(s)
heidel.diag(mc)

# Raftery-Lewis
# Want to the dependence factor less than 5.  If not the case,
# we'll need more samples
raftery.diag(mc)

# Geweke
# z scores for a test of equality (should be between -2 and 2)
# if converged
geweke.diag(mc)

