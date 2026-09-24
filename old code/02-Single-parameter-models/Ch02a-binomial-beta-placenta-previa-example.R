### Placenta previa example
# Placenta previa is a condition in which the placenta of an unborn child is
# implanted very low in the uterus, obstructing the child from
# a normal vaginal delivery. An early study concerning
# the sex of placenta previa births in Germany
# found that of a total of 980 births, 437 were female.
# How strong is the evidence that the proportion of female
# births in the population of placenta previa births is less than 0.485
# (the proportion in the general population)?

# Data distribution: y ~ Bin(n, theta) with n = 980 and observed y = 437
# Prior: theta ~ Beta(1, 1)
# Posterior: Beta(y + 1, n - y + 1) = Beta(438, 544)

#set parameters of Beta(a,b) posterior
a <- 438
b <- 544

mu <- a / (a + b) #determine mean of posterior
v <- (a * b) / (a + b)^2 / (a + b + 1) #determine variance of posterior

mu	#exact mean
qbeta(.5, a, b) #exact median
(ci0 <- qbeta(c(.025, .975), a, b)) #exact credible interval
pbeta(0.485, a, b) # exact posterior probability of interest

## Analysis using different conjugate prior distributions
# set data
# Data distribution: y ~ Bin(n, theta) with n = 980 and observed y = 437
# Prior: theta ~ Beta(alpha, beta)
# Posterior: Beta(y + alpha, n - y + beta)

# set data
n <- 980
y <- 437
# prior parameters
alpha <- c(1, 0.97, 2.425, 4.85, 9.7, 48.5, 97)
beta <- c(1, 1.03, 2.575, 5.15, 10.3, 51.5, 103)

#Determine corresponding posterior distribution parameters
a <-  alpha + y
b <- n - y + beta

a / (a + b) #exact means
qbeta(p = .5, a, b) #exact medians
cbind(qbeta(c(.025), a, b), qbeta(c(.975), a, b)) #exact posterior intervals
pbeta(0.485, a, b) # exact posterior probability of interest

# approximate posterior using simulation
x <- rbeta(1000, a, b) #1000 draws from posterior distribution
mean(x) #approximate mean of posterior
median(x) #approximate median of posterior
(ci1 <- quantile(x, prob = c(.025, .975))) #95% central posterior interval
mean(x <= 0.485) # approximate posterior probability

ci0
ci1

# plot results
theta <- seq(0, 1, by = 0.001) #create vector of possible thetas

#determine prior densities for each combination prior parameters
prior1 <- dbeta(theta, shape1 = alpha[1], shape2 = beta[1])
prior2 <- dbeta(theta, shape1 = alpha[2], shape2 = beta[2])
prior3 <- dbeta(theta, shape1 = alpha[3], shape2 = beta[3])
prior4 <- dbeta(theta, shape1 = alpha[4], shape2 = beta[4])
prior5 <- dbeta(theta, shape1 = alpha[5], shape2 = beta[5])
prior6 <- dbeta(theta, shape1 = alpha[6], shape2 = beta[6])
prior7 <- dbeta(theta, shape1 = alpha[7], shape2 = beta[7])

#determine posterior densities for each combination prior parameters
post1 <- dbeta(theta, shape1 = a[1], shape2 = b[1])
post2 <- dbeta(theta, shape1 = a[2], shape2 = b[2])
post3 <- dbeta(theta, shape1 = a[3], shape2 = b[3])
post4 <- dbeta(theta, shape1 = a[4], shape2 = b[4])
post5 <- dbeta(theta, shape1 = a[5], shape2 = b[5])
post6 <- dbeta(theta, shape1 = a[6], shape2 = b[6])
post7 <- dbeta(theta, shape1 = a[7], shape2 = b[7])

#construct 2x2 layout of plots
par(mfrow = c(2, 2))
#plot several prior/posterior combinations with legend and descriptive title
plot(theta, post1, type = "l", col = "orange")
lines(theta, prior1, col = "blue")
title(paste("alpha =", alpha[1], "beta =", beta[1]))
legend("topright", legend = c("posterior", "prior"),
       col = c("orange", "blue"), lwd = c(1, 1))

plot(theta, post2, type = "l", col = "orange")
lines(theta, prior2, col = "blue")
title(paste("alpha =", alpha[2], "beta =", beta[2]))
legend("topright", legend = c("posterior", "prior"),
       col = c("orange", "blue"), lwd = c(1, 1))

plot(theta, post5, type = "l", col = "orange")
lines(theta, prior5, col = "blue")
title(paste("alpha =", alpha[5], "beta =", beta[5]))
legend("topright", legend = c("posterior", "prior"),
       col = c("orange", "blue"), lwd = c(1, 1))

plot(theta, post7, type = "l", col = "orange")
lines(theta, prior7, col = "blue")
title(paste("alpha =", alpha[7], "beta =", beta[7]))
legend("topright", legend = c("posterior", "prior"),
       col = c("orange", "blue"), lwd = c(1, 1))
par(mfrow = c(1, 1))
