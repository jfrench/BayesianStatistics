### Example: Flexibility of the Beta prior

# Many conjugate priors are sufficiently flexible to approximately
# mimic your true prior beliefs.

# In this example, we look at the flexibility of the Beta
# distribution, which is a conjugate family for the
# Binomial family of Beta distributions.

# sequence of values at which to evaluate theta
theta = seq(0.001, 0.999, by = 0.001)

# Evaluate 3 beta densities with different shape parameters,
# (0.1, 0.1), (1, 1), and (50, 100)
d1 = dbeta(theta, 0.1, 0.1)
d2 = dbeta(theta, 1, 1)
d3 = dbeta(theta, 50, 100)

# plot the three different densities and note their distinct
# shapes
plot(theta, d3, type = "l", col = "orange", lty = 3,
     xlab = expression(theta),
     ylab = "density")
lines(theta, d2, col = "blue", lty = 1)
lines(theta, d1, col = "brown", lty = 2)
legend("topright",
       legend = c("Beta(1, 1)", "Beta(0.1, 0.1)", "Beta(50, 100)"),
       col = c("blue", "brown", "orange"),
       lwd = 1,
       lty = 1:3)


# a similar plot in ggplot2
# load ggplot2 if available and plot
if (requireNamespace("ggplot2")) {
  library(ggplot2)
  grp = rep(c("Beta(0.1, 0.1)", "Beta(1, 1)", "Beta(50, 100)"),
            each = length(theta))
  # convert data to data frame
  dtf = data.frame(theta = rep(theta, 3),
                   density = c(d1, d2, d3),
                   group = factor(grp))
  ggplot(data = dtf) +
    geom_line(aes(x = theta, y = density, col = group)) +
    theme_bw()
}

