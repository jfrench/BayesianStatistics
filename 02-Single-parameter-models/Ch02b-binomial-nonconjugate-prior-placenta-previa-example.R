###2.5 Placenta previa example
# Placenta previa is a condition in which the placenta of an unborn child is
# implanted very low in the uterus, obstructing the child from
# a normal vaginal delivery. An early study concerning
# the sex of placenta previa births in Germany
# found that of a total of 980 births, 437 were female.
# How strong is the evidence that the proportion of female
# births in the population of placenta previa births is less than 0.485
# (the proportion in the general population)?

# Data distribution: y ~ Bin(n, theta) with n = 980 and observed y = 437

# Analysis using a nonconjugate prior distribution

#create density function of prior
custom.prior <- function(x) {
  out <- numeric(length(x))

  for (i in 1:length(x)) {
    if (x[i] <=  .385) {
      out[i] <- .5
    } else if (x[i] >= .585) {
      out[i] <- .5
    } else if ((x[i] > .385) && (x[i] < .485)) {
      out[i] <- (.5 + 5 * (x[i] - .385)/.1)
    } else if ((x[i] >= .485) && (x[i] < .585)){
      out[i] <- (5.5 - 5 * (x[i] - .485)/.1)
    } else {
      out[i] <- 0
    }
  }
  return(out)
}

# theta values at which to evaluate densities
theta = seq(0, 1, by = 0.001)

#find custom prior density at specified values of theta
dcustom <- custom.prior(theta)

#plot nonconjugate density
par(mfrow = c(1, 1))
plot(theta, dcustom, type = "l",
     main = "nonconjugate prior density",
     xlab = "density")

#make sure prior integrates to 1
integrate(custom.prior, 0, 1, subdivisions = 1000)

#create function for posterior density
custom.post <- function(theta) {
  dbinom(x = 437, size = 980, prob = theta) *
    custom.prior(theta)/0.003616298
}

#make sure posterior integrates to 1
integrate(f = custom.post, 0, 1, subdivisions = 1000)

#plot posterior and prior distributions together with labels, legend
plot(theta, custom.post(theta), type = "l", col = "orange")
lines(theta, dcustom, col = "blue")
legend("topright", legend = c("posterior", "prior"),
       col = c("orange", "blue"), lwd = c(1, 1))

#create discrete approximation of posterior density (see p. 284)
discretep <- custom.post(theta)/sum(custom.post(theta))

plot(theta, discretep, col = "grey")

#create approximate posterior sample
theta.sample <- sample(theta, 100000, replace = TRUE, prob = discretep)

plot(density(theta.sample), main = "posterior simulations", xlab = "theta")
mean(theta.sample) #approximate mean
median(theta.sample) #approximate median
#approximate posterior interval
quantile(theta.sample, prob = c(.025, .975)) #approximate posterior interval

# comparing results
plot(theta, custom.post(theta), type = "l", col = "orange")
lines(theta, dcustom, col = "blue")
lines(density(theta.sample), col = "purple")
legend("topright", legend = c("posterior", "prior", "approximate posterior"),
       col = c("orange", "blue", "purple"), lwd = c(1, 1))



