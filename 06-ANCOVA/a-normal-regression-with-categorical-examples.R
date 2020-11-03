# generate some explanatory variables
x1 = sample(1:100, 20)
x2 = sample(1:100, 20)
# generate some data
a = 2 + 3 * x1 + rnorm(20, sd = 20)
b = 3 + 2 * x2 + rnorm(20, sd = 20)

# combine responses
y = c(a, b)
# combine explanatory variables
x = c(x1, x2)

# create group variable
group = rep(1:2, each = 20)

# constant mean model for all data (no groups)
par(mfrow = c(2, 3))
plot(x, y, pch = group, col = group + 2)
abline(h = mean(y))
title("Constant")

# simple linear regression model for all data (no groups)
plot(x, y, pch = group, col = group + 2)
abline(lm(y ~ x))
title("Single line")

# separate, constant mean for each group
plot(x, y, pch = group, col = group + 2)
abline(h = mean(a))
abline(h = mean(b))
title("Constant for each group")

# separate intercept for each group, common slope
lme = lm(y ~ x + I(rep(0:1, each = 20)))
plot(x, y, pch = group, col = group + 2)
abline(lme$coef[1], lme$coef[2])
abline(lme$coef[1] + lme$coef[3], lme$coef[2])
title("Parallel lines")

# common intercept, separate slope for each group
lmd = lm(y ~ x + I(rep(0:1, each = 20) * x))
plot(x, y, pch = group, col = group + 2)
abline(lmd$coef[1], lmd$coef[2])
abline(lmd$coef[1], lmd$coef[2] + lmd$coef[3])
title("Common intercept")

# separate intercept and slope for each group
lmf = lm(y ~ x + rep(0:1, each = 20) + I(rep(0:1, each = 20) * x))
plot(x, y, pch = group, col = group + 2)
abline(lmf$coef[1], lmf$coef[2])
abline(lmf$coef[1]  + lmf$coef[3] , lmf$coef[2] + lmf$coef[4])
title("Separate lines")