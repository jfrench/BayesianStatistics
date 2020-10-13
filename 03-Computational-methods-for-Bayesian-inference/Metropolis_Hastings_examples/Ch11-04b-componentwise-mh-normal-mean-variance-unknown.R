### Create componentwise metropolis hastings
### sampler for normal(mu, sigma^2) sampling distribution
### with mu and sigma^2 unknown.

# Data distribution: y1, y2, ..., yn ~ iid N(mu, sigma^2).
# Prior distribution: 
# mu | sigma^2 ~ N(mu0, sigma^2/k0) w/ mu0 = 1.9 and k0 = 1
# sigma^2 ~ Inv-Chisquare(nu0, sigma0^2) with nu0 = 1 and sigma0 = 0.1.

# Jumping distributions:
# mu.star | mu^(t-1) ~ N(mu^(t-1), mustar.sigma^2) w/ mustar.sigma^2 = 0.2^2
# sigmasq.star ~ Inv-chisq(df.sigmasqstar, scale.sigmasqstar^2)) with 
# df.sigmasqstar = 2 and scale.sigmasqstar = 0.1.

### Midge Example
### Chapter 5 of A First Course in Bayesian Statistical Methods by PD Hoff

#Grogan and Wirth (1981) provide data on the wing length in millimeters 
#of nine members of a species of midge (small, two-winged flies). 
#From these nine measurements we wish to make inference on the population mean
#theta. 

#Studies of other populations suggest that the true mean should be around 1.9 mm
#with a standard deviation of 0.1.  However, this population may be different
#from the others, so we choose k0 and vu0 = 1 so that the prior distributions
#are only weakly centered around these estimates from other populations.

# load useful distributions
source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")
library(coda)

#Set prior parameters
mu0 <- 1.9
k0 <- 1
nu0 <- 1
sigma0 <- 0.1

#define number of simulations, data, sample size, sample st dev, sample mean
B <- 100000
y <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n <- length(y)

#Initial value for sigma
s <- sd(y)
#Initial value for mu
m <- mean(y)

cmh = function(B, start, jump.parm1, jump.parm2) {
  mustar.sigma <- jump.parm1
  df.sigmasqstar <- jump.parm2[1]
  scale.sigmasqstar <- jump.parm2[2]
  
  theta = matrix(0, nrow = B + 1, ncol = length(start))
  theta[1, ] = start

  for (i in 2:(B + 1)) {
    # component for mu
    mu.star = rnorm(1, mean = theta[i - 1, 1], sd = mustar.sigma)
    # print(paste("mu", mu.star)) used for debugging
    # what are the current values during a crash?
    
    # notice where we use the updated mu and where we do not
    # the sigmasq value just uses the previous iteration sigmasq
    num.logr <- sum(log(dnorm(y, mean = mu.star, sd = sqrt(theta[i - 1, 2])))) + 
      log(dnorm(mu.star, mean = mu0, sd = sqrt(theta[i - 1, 2]/k0))) + 
      log(dinvchisq(theta[i - 1, 2], df = nu0, scale = sigma0)) - 
      log(dnorm(mu.star, mean = theta[i - 1, 1], sd = mustar.sigma)) 
    den.logr <- sum(log(dnorm(y, mean = theta[i - 1, 1], sd = sqrt(theta[i - 1, 2])))) + 
      log(dnorm(theta[i - 1, 1], mean = mu0, sd = sqrt(theta[i - 1, 2]/k0))) + 
      log(dinvchisq(theta[i - 1, 2], df = nu0, scale = sigma0)) - 
      log(dnorm(theta[i - 1, 1], mean = mu.star, sd = mustar.sigma)) 
      
    logr <- num.logr - den.logr
    
    if (log(runif(1)) <= min(logr, 0)) {
      theta[i, 1] <- mu.star
    } else {
      theta[i, 1] <- theta[i - 1, 1]
    }

    # component for sigmasq
    # use mu for current iteration
    # change value of sigmasq in numerator and denominator
    
    sigmasq.star = rinvchisq(1, df = df.sigmasqstar, scale = scale.sigmasqstar)
    # print(sigmasq.star) used for debugging
    # what are the current values during a crash?
    
    # notice where we use the updated mu and where we do not
    # the sigmasq value just uses the previous iteration sigmasq
    num.logr <- sum(log(dnorm(y, mean = theta[i, 1], sd = sqrt(sigmasq.star)))) + 
      log(dnorm(theta[i, 1], mean = mu0, sd = sqrt(sigmasq.star/k0))) + 
      log(dinvchisq(sigmasq.star, df = nu0, scale = sigma0)) - 
      log(dinvchisq(sigmasq.star, df = df.sigmasqstar, scale = scale.sigmasqstar)) 
    den.logr <- sum(log(dnorm(y, mean = theta[i, 1], sd = sqrt(theta[i - 1, 2])))) + 
      log(dnorm(theta[i, 1], mean = mu0, sd = sqrt(theta[i - 1, 2]/k0))) + 
      log(dinvchisq(theta[i - 1, 2], df = nu0, scale = sigma0)) - 
      log(dinvchisq(theta[i - 1, 2], df = df.sigmasqstar, scale = scale.sigmasqstar)) 
    logr <- num.logr - den.logr

    if (log(runif(1)) <= min(logr, 0)) {
      theta[i, 2] <- sigmasq.star
    } else {
      theta[i, 2] <- theta[i - 1, 2]
    }
  }
  return(theta)
}

# run 3 chains with different starting values
# the sigmasq parameter doesn't mix very well.
B = 100000 # run with 10000 and 100000
# mu.sigma = .2
# df.sigmasq = 2
# scale.sigmasq = .1
chain1 = cmh(B, start = c(m, s^2), jump.parm1 = .2, jump.parm2 = c(2, .1))
chain2 = cmh(B, start = c(1, 0.5^2), jump.parm1 = .2, jump.parm2 = c(2, .1))
chain3 = cmh(B, start = c(3, 0.05^2), jump.parm1 = .2, jump.parm2 = c(2, .1))

# values to keep
keep = seq(B/2 + 1, B + 1, by = 1)
# combine chains into mcmc list.  Discard half as warmup. 
mc <- mcmc.list(mcmc(chain1[keep, ]),
                mcmc(chain2[keep, ]),
                mcmc(chain3[keep, ]))

# trace plot for the 3 chains
traceplot(mc)

# acf plot
autocorr.plot(mc)

# effective sample size
effectiveSize(mc)

# Calculate scale reduction factor
gelman.diag(mc, autoburnin = FALSE)

# Plot Gelman statistic
gelman.plot(mc, autoburnin = FALSE)

# Heidelberg & Welch
heidel.diag(mc) # if the halfwidth test fails, extend the chain(s)

# Raftery-Lewis
# like to see dependence factor less than 6.  If not the case,
# we'll need more samples
raftery.diag(mc) 

# Geweke
geweke.diag(mc) # z scores for a test of equality (should be between -2 and 2)
# if converged


