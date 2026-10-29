
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] d2;
  real<lower=0> v;
}
parameters {
  real beta0;
  real beta1;
  real alpha2;
  real delta2;
  real<lower=0> sigmasq;
}
transformed parameters {
  vector[n] mu = beta0 + beta1*x + alpha2*d2 + delta2*(x .* d2);
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  delta2 ~ normal(0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);
  y ~ normal(mu, sqrt(sigmasq));
}
generated quantities {
  real<lower=0,upper=1> R2 = variance(mu)/(variance(mu)+sigmasq);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sqrt(sigmasq));
}

