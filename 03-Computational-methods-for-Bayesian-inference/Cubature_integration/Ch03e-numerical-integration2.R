 # Log-odds example
# Binomial-Beta sampling/prior distribution combo for proportion theta
# Posterior distribution for theta is Beta(442, 420)
# Find posterior distribution for gamma = log(theta / (1 - theta))

# Specify number of samples to take
B <- 10000

# Sample nsim samples of theta from its posterior distribution
theta <- rbeta(B, 442, 420)

# Apply the logit transformation to posterior theta samples to get samples from gamma
# posterior distribution
gamma <- log(theta / (1 - theta))

# Estimated mean of gamma posterior
mean(gamma)
# Estimated variance of gamma posterior
var(gamma)
# Estimated standard error of sample mean of gamma posterior
sqrt(var(gamma)/B)
# 95% posterior interval for gamma
quantile(gamma, prob = c(.025, .975))

# Empirical posterior distribution of gamma
# Apply a kernel density smoother to sample to obtain empirical density
plot(density(gamma), main = expression(paste("Approximate posterior of ", gamma)),
     xlab = expression(gamma))