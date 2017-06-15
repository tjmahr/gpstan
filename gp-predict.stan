// Predict from Gaussian Process
// All data parameters must be passed as a list to the Stan call
// Based on original file from https://code.google.com/p/stan/source/browse/src/models/misc/gaussian-process/

data {
  int<lower=1> N1;
  vector[N1] x1;
  vector[N1] y1;
  int<lower=1> N2;
  vector[N2] x2;
  real sigma_sq;
  real eta_sq;
  real rho_sq;
}
transformed data {
  int<lower=1> N;
  vector[N1 + N2] x;
  vector[N1 + N2] mu;
  cov_matrix[N1 + N2] Sigma;
  N = N1 + N2;
  for (n in 1:N1) x[n] = x1[n];
  for (n in 1:N2) x[N1 + n] = x2[n];
  for (i in 1:N) mu[i] = 0;

  for (i in 1:N)
    for (j in 1:N)
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j], 2))
                    + i == j ? sigma_sq : 0.0;
}
parameters {
  vector[N2] y2;
}
model {
  vector[N] y;
  for (n in 1:N1) y[n] = y1[n];
  for (n in 1:N2) y[N1 + n] = y2[n];

  y ~ multi_normal(mu, Sigma);
}

