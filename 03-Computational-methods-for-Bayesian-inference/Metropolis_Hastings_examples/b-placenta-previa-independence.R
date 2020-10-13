library(coda)

### Placenta previa example

# construction Metropolis-Hastings algorithm with an independence
# sampler to to estimate the posterior of the proportion of female
# births among placenta previa patients

# Data distribution: y | theta ~ Binomial(n, theta)
# Prior: theta ~ U(0, 1)
# n = 980, observed y = 437

# Proposal distribution: theta.star ~ Beta(alpha, beta)

mh = function(B, start, jump.parm) {
  alpha = jump.parm[1]
  beta = jump.parm[2]

  theta = numeric(B + 1)
  theta[1] = start

  for (i in 2:length(theta)) {
		theta.star = rbeta(1, alpha, beta)
		num.logr = log(dbinom(437, size = 980, prob = theta.star)) +
             log(dbeta(theta.star, 1, 1)) -
             log(dbeta(theta.star, alpha, beta))
		den.logr = log(dbinom(437, size = 980, prob = theta[i - 1])) +
             log(dbeta(theta[i - 1], 1, 1)) -
             log(dbeta(theta[i - 1], alpha, beta))
		logr = num.logr - den.logr

		if (log(runif(1)) <= min(logr, 0)) {
			theta[i] = theta.star
		} else {
			theta[i] = theta[i - 1]
		}
	}
	return(theta)
}

set.seed(77)
# length of chain
B = 10000
# values to keep
keep = (B/2 + 1):(B + 1)
# jump parameters.  Compare for c(1, 1), c(2, 2), c(3, 3)
jpar = c(1, 1)
chain1 = mh(B, start = 0.01, jump.parm = jpar)
chain2 = mh(B, start = 0.25, jump.parm = jpar)
chain3 = mh(B, start = 0.5, jump.parm = jpar)
chain4 = mh(B, start = 0.75, jump.parm = jpar)
# chain5 = mh(B, start = 0.99, jump.parm = jpar) # bad starting value
chain5 = mh(B, start = 0.9, jump.parm = jpar)

library(coda)

# combine chains into mcmc list.  Discard half as warmup
mc = mcmc.list(mcmc(chain1[keep]),
                mcmc(chain2[keep]),
                mcmc(chain3[keep]),
                mcmc(chain4[keep]),
                mcmc(chain5[keep]))

# trace plot for the 5 chains
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

### An example of thinning (keep every 10th observation)
# length of chain
B = 10000*10
# values to keep
keep = seq(B/2 + 1, B + 1, by = 10)
# jump parameters.  Compare for c(1, 1), c(2, 2), c(3, 3)
jpar = c(3, 3)
chain1 = mh(B, start = 0.01, jump.parm = jpar)
chain2 = mh(B, start = 0.25, jump.parm = jpar)
chain3 = mh(B, start = 0.5, jump.parm = jpar)
chain4 = mh(B, start = 0.75, jump.parm = jpar)
chain5 = mh(B, start = 0.95, jump.parm = jpar)

# combine chains into mcmc list.  Discard half as warmup.  Keep every
# 10th observation
mc = mcmc.list(mcmc(chain1[keep]),
                mcmc(chain2[keep]),
                mcmc(chain3[keep]),
                mcmc(chain4[keep]),
                mcmc(chain5[keep]))

# trace plot for the 5 chains
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
