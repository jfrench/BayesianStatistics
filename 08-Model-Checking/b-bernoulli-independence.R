# use bayesplot for posterior predictive check
library(bayesplot)

# Check independence of binomial trials

# Data distribution: y|theta ~ Binomial(n, theta)
# Prior distribution: theta ~ U(0, 1)
# Postieor distribution:
# theta |y ~ Beta(sum(y) + 1, n - sum(y) + 1)

#observed data
y = c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

# other information
sumy = sum(y)
n = length(y)
B = 1e4 #10000 samples

# sample from posterior distribution of theta
theta = rbeta(B, sumy + 1, n - sumy + 1)

# number of switches between 0 and 1 using rle
# (run length encoding) function
nswitch = function(x) length(rle(x)$values) - 1

# generate yrep
yrep = matrix(0, nrow = B, ncol = n)
for (i in 1:B) {
  # sample n bernoulli trials with success probability theta[i]
  yrep[i, ] = rbinom(n, size = 1, prob = theta[i])
}

nswitch_y = nswitch(y)
nswitch_yrep = apply(yrep, 1, nswitch)

# posterior predictive probability
(sum(nswitch_yrep >= nswitch_y) + 1)/(B + 1)

# histogram of number of switches for yrep
# compared to observed value
hist(nswitch_yrep, xlab = "switches", main = "")
abline(v = nswitch_y)

ppc_stat(y, yrep, stat = nswitch)
