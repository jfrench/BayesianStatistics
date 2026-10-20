
data {
  int<lower=1> n; // number of births
  int<lower=0> y; // number of female placenta previa births
  real<lower=0> alpha; // prior parameters
  real<lower=0> beta; // prior parameters
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  y ~ binomial(n, theta);  //data distribution
  theta ~ beta(alpha, beta); //prior distribution
}
generated quantities{
  int<lower=0> ytilde;
  ytilde = binomial_rng(n, theta);
}

