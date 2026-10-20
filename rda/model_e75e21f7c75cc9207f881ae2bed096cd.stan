
data {
  int<lower=1> n; // number of observations
  vector[n] y; // measurements
}
parameters {
  real logsigma;
  real mu;
}
transformed parameters {
  real<lower=0> sigma = exp(logsigma);
}
model {
  y ~ normal(mu, sigma);
  // Flat improper priors for mu and logsigma
}
generated quantities {
  real<lower=0> sigmasq = square(sigma);
  vector[n] yrep;
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sigma);
  }
}

