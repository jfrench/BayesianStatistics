
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] tutor2;
  vector[n] tutor3;
  vector[n] tutor4;
}
parameters {
  real<lower=0> sigmasq;
  real mu0;
  real alpha2;
  real alpha3;
  real alpha4;
}
model {
  vector[n] mu = mu0 + alpha2*tutor2 + alpha3*tutor3 + alpha4*tutor4;
  mu0 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  alpha3 ~ normal(0, 100);
  alpha4 ~ normal(0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);
  y ~ normal(mu, sqrt(sigmasq));
}
generated quantities {
  vector[4] group_mean;
  group_mean[1] = mu0;
  group_mean[2] = mu0 + alpha2;
  group_mean[3] = mu0 + alpha3;
  group_mean[4] = mu0 + alpha4;
}

