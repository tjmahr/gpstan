// Sample from Gaussian process
// All data parameters must be passed as a list to the Stan call
// Based on original file from https://code.google.com/p/stan/source/browse/src/models/misc/gaussian-process/

data {
  int<lower=1> N;
  real x[N];
  real eta_sq;
  real rho_sq;
  real sigma_sq;
}
transformed data {
  vector[N] mu;
  cov_matrix[N] Sigma;
  for (i in 1:N)
    mu[i] = 0;
  for (i in 1:N)
    for (j in 1:N)
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j], 2)) + i == j ? sigma_sq : 0.0;
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu,Sigma);
}
