# Summarizing posterior inference

# load and install TeachingDemos package if necessary
if (!require(TeachingDemos)) {
  install.packages("TeachingDemos",
                   repos = "https://cran.r-project.org")
  library(TeachingDemos)
}

#posterior distribution for theta is Beta(10, 2)
#set shape parameters
a = 10
b = 2

#create sequence of possible values of theta
theta = seq(0, 1, len = 1000)
#find posterior density for these values of theta
d = dbeta(theta, shape1 = a, shape2 = b)

# plot of density
plot(theta, d, type = "l",
     ylab = "density",
     xlab = expression(theta))

#find 95% central credible interval and hpd
cci = qbeta(c(.025, .975), shape1 = a, shape2 = b)
cci
hpdi = hpd(qbeta, shape1 = a, shape2 = b, conf = .95)
hpdi

#draw vertical lines on plot for cpi and hpdi.  Add legend
abline(v = cci, col = "blue")
abline(v = hpdi, col = "orange")
legend("topleft", legend = c("central", "HPD"), lty = 1,
       col = c("blue", "orange"))

# MAP estimate
optimize(f = dbeta, interval = 0:1, shape1 = 10, shape2 = 2,
         maximum = TRUE)



