# An example of finding the posterior mode
# (maximum a posteriori estimation)

# Data distribution: y_1, ..., y_n | \theta ~ i.i.d Poisson(theta)
# Prior distribution: theta ~ Exp(1)

# The log posterior is proportional to
# const -(n + 1) * theta + sum(y) * log(theta)

# Generate some synthetic data
# poisson with mean of 0.5
# in practice, we just have some data
set.seed(8)
y = rpois(10, lambda = 0.5)

# create function for log unnormalized posterior
# theta is unknown parameter
# sumy = sum(y)
# n = sample size
lup = function(theta, sumy, n) {
  -(n + 1) * theta + sumy * log(theta)
}

# use optimize function to optimize lup between 0
# and 1. Start interval just above 0 to avoid
# numerical issues. Set maximum to TRUE.
# sumy and n must be passed as arguments
(map = optimize(f = lup, interval = c(0.0001, 10), maximum = TRUE,
                sumy = sum(y), n = length(y)))
# maximum occurs at theta = 0.3636

# plot results to double-check accuracy
# evaluate lup at many values of theta
mytheta = seq(0.0001, 1, length = 1000)
# the Vectorize function "vectorizes" a function
# that normally only takes a single value for
# an argument.
# We will Vectorize lup so that we can evaluate
# a sequence of theta values
vlup = Vectorize(lup, vectorize.args = "theta")
eval_lup = vlup(mytheta, sumy = sum(y), n = length(y))

# plot function
# show maximum
plot(mytheta, eval_lup, type = "l",
     xlab = "theta", ylab = "propto log posterior")
abline(v = map$maximum)

# A second approach without the simplification
# since log (p(y|theta) * p(theta))
# = log(p(y|theta)) + log(p(theta))
# = log(prod(p(yi | theta))) + log(p(theta))
# = sum(log(p(yi | theta))) + log(p(theta))
lup2 = function(theta, y) {
  sum(dpois(y, lambda = theta, log = TRUE)) +
    dexp(theta, rate = 1, log = TRUE)
}

# same results!
optimize(f = lup2, interval = c(0.0001, 10), maximum = TRUE,
         y = y)

