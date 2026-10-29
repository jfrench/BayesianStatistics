
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] cases;
  vector[n] dist;
  real prior_mean;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_beta1_sd;
  real<lower=0> prior_beta2_sd;
  real<lower=0> prior_sigma_rate;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}
model {
  beta0 ~ normal(prior_mean, prior_intercept_sd);
  beta1 ~ normal(0, prior_beta1_sd);
  beta2 ~ normal(0, prior_beta2_sd);
  sigma ~ exponential(prior_sigma_rate);
  y ~ normal(beta0 + beta1*cases + beta2*dist, sigma);
}
generated quantities {
  real sigmasq = square(sigma);
  vector[n] mu = beta0 + beta1*cases + beta2*dist;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
}

