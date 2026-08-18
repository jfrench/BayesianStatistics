# Example:  Example 1 of Chapter 9 of Bayesian Modeling
# Using Winbugs by Ntzoufras
#
# Consider a problem related to blood pressure measurements.  We consider
# the repeated measurements of blood pressure from 20 health individuals.
# Our aim is to estimate both the between-subject and within-subject variability.
#
# We assume that is y_ij = mu + a_i + epsilon_ij
# We consider a hierarchical model
# a_i ~ N(0, sigmasq_a), epsilon_ij ~ N(0, sigmasq)
#
# Equivalently, Y_ij ~ N(mu_i, sigma^2) with mui = mu + a_i
# and a_i ~ N(0, sigmasq_a)
#
# We assume low information prior distributions for the regression coefficients,
# mu ~ N(0, 1000), sigmasq ~ IG(0.001, 0.001), sigmasq_a ~ IG(0.001, 0.001)

library(brms)

# Enter data manually
n <- 20
y <- c(108, 98, 91, 94, 93, 96, 104, 99, 99, 97, 95,
             98, 93, 97, 99, 96, 90, 100, 92, 95, 101, 89,
             97, 97, 97, 100, 96, 95, 106, 100, 100, 98,
             90, 99, 88, 98, 92, 92, 100, 101)
bdat <- data.frame(y = y,
                   subject = factor(rep(1:20, each = 2)),
                   measurement = rep(c(1, 2), times = 20))
bmod <- brm(y ~ 1 + (1 | subject),
            data = bdat,
            family = gaussian(),
            prior = c(set_prior("normal(0, sqrt(1000))", class = "Intercept"),
                      set_prior("student_t(4, 0, 1)", class = "sd"),
                      set_prior("student_t(4, 0, 1)", class = "sd", group = "subject")),
            iter = 1)
bsamples_brms <- update(bmod,
                        chains = 2,
                        iter = 1e+4, seed = 23)


bmod2 <- stan_glmer(y ~ 1 + (1 | subject),
            data = bdat,
            family = gaussian(),
            prior_intercept = normal(0, sqrt(1000),
                                     autoscale = FALSE),
            prior_aux = student_t(4, 0, 1, autoscale = FALSE),
            chains = 4, iter = 5e+4, seed = 23)
summary(bmod2, pars = c("(Intercept)", "sigma",
                        "Sigma[subject:(Intercept),(Intercept)]"))


# check convergence
sums_bmod <- summary(bsamples_brms)

$summary[,"Rhat"]

# posterior summary
post_sum = summary(blood_samples, probs = c(0.025, 0.975))$summary
round(post_sum[c("mu", "sigmasq", "sigmasq_a", "tvar", "cor",
                 "sigma", "sigma_a"), c("mean", "sd", "2.5%", "97.5%")],
      2)
