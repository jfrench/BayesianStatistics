
data {
  int<lower=1> n;
  vector[n] y;
  real mu0;
  real<lower=0> k0;
  real<lower=0> nu0;
  real<lower=0> sigma0;
}
parameters {
  real<lower=0> sigmasq;
  real mu;
}
model {
  y ~ normal(mu, sqrt(sigmasq));
  mu ~ normal(mu0, sqrt(sigmasq / k0));
  sigmasq ~ inv_gamma(nu0 / 2, nu0 * square(sigma0) / 2);
}
generated quantities {
  vector[n] yrep;
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sqrt(sigmasq));
  }
}

