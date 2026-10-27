
data {
  int<lower=1> n;
  vector[n] y;
  matrix[n,3] X;
  vector[3] mu0;
  cov_matrix[3] V;
  real<lower=0> g;
  real<lower=0> v;
}
parameters {
  vector[3] beta;
  real<lower=0> prec;
}
transformed parameters {
  real<lower=0> sigmasq = 1 / prec;
}
model {
  beta ~ multi_normal(mu0, sigmasq * g * V);
  prec ~ gamma(0.01, 0.01);
  y ~ normal(X * beta, sqrt(sigmasq));
}
generated quantities {
  real Rbsq = 1 - sigmasq / v;
  vector[n] mu = X * beta;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
}

