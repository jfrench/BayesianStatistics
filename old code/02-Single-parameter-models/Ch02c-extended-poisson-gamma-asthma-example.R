### Poisson-gamma asthma example
### Taken from Chapter 2 of Bayesian Data Analysis,
### 3rd edition, by Gelman et al.

# Out of a population of 200,000 people,
# 3 persons died of asthma last year.
# A crude estimate of asthma mortality rate
# in the city is 1.5 cases per 100,000 persons
# per year. The data distribution of y, the
# number of deaths in a city of 200,000 in
# one year, may be expressed as Poisson(2.0 * theta),
# where theta represents the true underlying
# mortality rate in the city (in cases per
# 100,000 persons per year). Our observation y=3
# with an exposure rate x=2.0.

# data distribution: y|theta ~ Poisson(2 * theta)
# prior: theta ~ Gamma(3, 5)
# posterior: theta|y ~ Gamma(1 + )

#setup prior
alpha <- 3
beta <- 5

#verify that a Gamma(3, 5) random variable is almost always
#less than 1.5
pgamma(1.5, shape = alpha, rate = beta)

# 100000 draws from posterior
post.draws <- rgamma(100000, shape = alpha + 3, rate = beta + 2)
mean(post.draws) #approximate mean

#draw, label, and title histogram of draws
plot(density(post.draws), xlab = "theta", main = "100000 draws from Gamma(6,7)")
theta = seq(0.0001, 3, len = 1000)
post1 = dgamma(theta, shape = alpha + 3, rate = beta + 2)
plot(theta, post1,
     main = "Gamma(6, 7) density", type = "l",
     ylab = "density")

#posterior with more data

#100000 draws from posterior
post.draws2 <- rgamma(100000, shape = alpha + 30, rate = beta + 2 * 10)

# Example continued
# Suppose that we have ten years of data. We observe 30 deaths
# over the 10 years. Assuming a constant population for each
# year, determine the posterior

#draw, label, and title density of draws
plot(density(post.draws2), xlab = "theta", main = "100000 draws from Gamma(33,25)")
mean(post.draws2)
1 - pgamma(1, shape = alpha + 30, rate = beta + 2 * 10)
theta = seq(0.0001, 3, len = 1000)
post2 = dgamma(theta, shape = alpha + 30, rate = beta + 2 * 10)
plot(theta, post2,
     main = "Gamma(33, 25) density", type = "l",
     ylab = "density",
     xlim = c(0.5, 2.5))

# overlay densities side by side
posterior <- factor(rep(c("Gamma(6, 7)", "Gamma(33, 25)"), each = length(theta)))
df <- data.frame(theta = rep(theta, 2),
                 density = c(post1, post2),
                 posterior = posterior)
library(ggplot2)
ggplot(df, aes(x = theta, y = density, group = posterior)) +
  geom_path(aes(col = posterior)) + theme_bw()

library(lattice)
xyplot(density ~ theta, data = df, group = posterior,
       type = "l",
       auto.key = TRUE)
