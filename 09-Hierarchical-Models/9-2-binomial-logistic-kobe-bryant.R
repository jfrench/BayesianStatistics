# Example 9.2 from Bayesian Modeling Using WinBugs
# Consider the probability of field goal success for Kobe Bryant
# A simple hierarchical model
# Data distribution: y_i ~ Binomial(N_i, pi_i)
# Prior distributions:
# logit(pi_i) = mu_i
# mu_i ~ i.i.d., N(mu_theta, sigmasq_theta)
# mu_theta ~ N(0, 1000)
# sigmasq_theta ~ Inv-gamma(0.01, 0.01)
# A state space model assumes:
# mu_i ~ N(mu_{i-1}, sigmasq_theta)

kobe <-
structure(list(SEASON = structure(as.integer(c(4, 5, 6, 7, 8,
9, 10, 11)), .Label = c(" 1996-97", " 1997-98", " 1998-99", " 1999-00",
" 2000-01", " 2001-02", " 2002-03", " 2003-04", " 2004-05", " 2005-06",
" 2006-07"), class = "factor"), TEAM = structure(as.integer(c(1,
1, 1, 1, 1, 1, 1, 1)), .Label = "LAL     ", class = "factor"),
    GAMES = c(66, 68, 80, 82, 65, 66, 80, 42), FGTOTAL = c(554,
    701, 749, 868, 516, 573, 978, 399), FGTOTATT = c(1183, 1510,
    1597, 1924, 1178, 1324, 2173, 845)), .Names = c("SEASON",
"TEAM", "GAMES", "FGTOTAL", "FGTOTATT"), row.names = c("4", "5",
"6", "7", "8", "9", "10", "11"), class = "data.frame")

library(rstan)

# fixed effects model
femod = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // field goals made
  int N[n];  // field goals attempted
}
parameters {
  vector[n] theta;  // probability of success
}
model {
  // prior distributions
  for(i in 1:n) theta[i] ~ normal(0, sqrt(1000));
  // data distribution
  y ~ binomial_logit(N, theta);
}
generated quantities {
  real mu_theta;
  real sigma_theta;
  real pi_mu_theta;
  real pi[n];
  vector[n] log_lik;
  mu_theta = mean(theta);
  sigma_theta = sd(theta);
  pi_mu_theta= 1/(1 + exp(-mu_theta))*100;
  for(i in 1:n) {
    pi[i] = 1/(1 + exp(-theta[i]))*100;
    log_lik[i] = binomial_lcdf(y[i] | N[i], pi[i]/100);
  }
}
"

# simple hierarchical model
shmod = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // field goals made
  int N[n];  // field goals attempted
}
parameters {
  real mu_theta;          // common mean
  real<lower=0> sigmasq_theta;  // variance of within-subject random effects
  vector[n] theta;  // probability of success
}
model {
  // prior distributions
  mu_theta ~ normal(0, 10);
  sigmasq_theta ~ inv_gamma(0.01, 0.01);
  for(i in 1:n) theta[i] ~ normal(mu_theta, sqrt(sigmasq_theta));

  // data distribution
  y ~ binomial_logit(N, theta);
}
generated quantities {
  real pi_mu_theta;
  real sigma_theta;
  real pi[n];
  vector[n] log_lik;
  sigma_theta = sqrt(sigmasq_theta);
  pi_mu_theta = 1/(1 + exp(-mu_theta))*100;
  for(i in 1:n) {
    pi[i] = 1/(1 + exp(-theta[i]))*100;
    log_lik[i] = binomial_lcdf(y[i] | N[i], pi[i]/100);
  }
}
"

# state space model
ssmod = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // field goals made
  int N[n];  // field goals attempted
}
parameters {
  real<lower=0> sigmasq_theta;  // variance of within-subject random effects
  vector[n] theta;  // probability of success
}
model {
  // prior distributions
  sigmasq_theta ~ inv_gamma(0.01, 0.01);
  theta[1] ~ normal(0, 10);
  for(i in 2:n) theta[i] ~ normal(theta[i-1], sqrt(sigmasq_theta));

  // data distribution
  y ~ binomial_logit(N, theta);
}
generated quantities {
  real pi[n];
  vector[n] log_lik;
  for(i in 1:n) {
    pi[i] = 1/(1 + exp(-theta[i]))*100;
    log_lik[i] = binomial_lcdf(y[i] | N[i], pi[i]/100);
  }
}
"

# create the data list
kobe_data = list(n = nrow(kobe), y = kobe$FGTOTAL, N = kobe$FGTOTATT)

if(!file.exists("example_9_2.rda")) {
# draw samples from the models
fe_mod <- stan(model_code = femod, data = kobe_data,
               iter = 5e4, seed = 32)
sh_mod <- stan(model_code = shmod, data = kobe_data,
                      iter = 5e4, seed = 24)
ss_mod <- stan(model_code = ssmod, data = kobe_data,
               iter = 5e4, seed = 8)
save(fe_mod, sh_mod, ss_mod, file = "example_9_2.rda")
}
load(file = "example_9_2.rda")

# posterior summaries
round(summary(fe_mod)$summary[c(9:19), c("mean", "sd", "2.5%", "97.5%")], digits = 1)
round(summary(sh_mod)$summary[c(1, 11:20), c("mean", "sd", "2.5%", "97.5%")], digits = 1)
round(summary(ss_mod)$summary[c(10:17), c("mean", "sd", "2.5%", "97.5%")], digits = 1)

# display 95% posterior interval for pi
# for each model
# get median, 95% posterior interval each each pi[i]
df1 = as.data.frame(round(summary(fe_mod)$summary[c(12:19), c("50%", "2.5%", "97.5%")], digits = 1))
df2 = as.data.frame(round(summary(sh_mod)$summary[c(13:20), c("50%", "2.5%", "97.5%")], digits = 1))
df3 = as.data.frame(round(summary(ss_mod)$summary[c(10:17), c("50%", "2.5%", "97.5%")], digits = 1))
# combine information
df = as.data.frame(rbind(df1, df2, df3))
names(df) = c("median", "lower", "upper")
# indicate model
df$model = rep(c("FE", "SH", "SS"), each = 8)
# decide x-axis position of each interval
sq = 1:8
df$x = c(sq, sq - 0.25, sq + 0.25)
# for naming intervals
s = rep(c("1999-00",
          "2000-01",
          "2001-02",
          "2002-03",
          "2003-04",
          "2004-05",
          "2005-06",
          "2006-07"), each = 1)

# plot results
ggplot(df, aes(x = x, lty = model, col = model)) +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  geom_point(aes(x = x, y = median)) +
  ggtitle("posterior summaries of field goal success") +
  ylab("field goal success percentage") +
  xlab("season") +
  scale_x_continuous(breaks = 1:8, labels = s) +
  theme(axis.text.x = element_text(angle = 90))

library(loo)
# extract log likelihoods
ll_fe = extract_log_lik(fe_mod, merge_chains = FALSE)
ll_sh = extract_log_lik(sh_mod, merge_chains = FALSE)
ll_ss = extract_log_lik(ss_mod, merge_chains = FALSE)


# compute waic on fitted models
waic_fe = waic(ll_fe)
waic_sh = waic(ll_sh)
waic_ss = waic(ll_ss)

# compute relative efficiency of log likelihoods
r_eff_fe = exp(relative_eff(ll_fe))
r_eff_sh = exp(relative_eff(ll_sh))
r_eff_ss = exp(relative_eff(ll_ss))

# compute looic for each model
looic_fe = loo(ll_fe, r_eff = r_eff_fe)
looic_sh = loo(ll_sh, r_eff = r_eff_sh)
looic_ss = loo(ll_ss, r_eff = r_eff_ss)

# compare results for three models
compare(waic_fe, waic_sh, waic_ss)
compare(loo_fe, loo_sh, loo_ss)


