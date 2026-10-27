
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] cases;
  vector[n] dist;
  real<lower=0> v;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> prec;
}
transformed parameters {
  real<lower=0> sigma = inv_sqrt(prec);
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  prec ~ gamma(0.01, 0.01);
  y ~ normal(beta0 + beta1 * cases + beta2 * dist, sigma);
}
generated quantities {
  real<lower=0> sigmasq = 1 / prec;
  real Rbsq = 1 - sigmasq / v;
  vector[n] mu = beta0 + beta1 * cases + beta2 * dist;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
}

