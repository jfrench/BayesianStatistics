
data {
  int<lower=1> n; // number of trials
  array[n] vector[2] y; // 2d array of size n
  vector[2] mu0; // prior mean for mu
  real<lower=0> k0; // prior number of observations for mu
  real<lower=0> nu0; // df for prior for Sigma
  cov_matrix[2] L0; // scale of prior for Sigma
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
}
model {
  // specify my prior distributions for mu and Sigma
  mu ~ multi_normal(mu0, Sigma/k0);
  Sigma ~ inv_wishart(nu0, L0);
  // specify the data distribution for each observation
  for(i in 1:n){
    y[i] ~ multi_normal(mu, Sigma);
  }
}
generated quantities {
  vector[2] ytilde; //pre- and post-test measurements for new students
  ytilde = multi_normal_rng(mu, Sigma);
}

