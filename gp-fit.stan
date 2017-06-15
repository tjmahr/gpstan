# Estimate from Gaussian Process from data
data {
  # x1 and y1 are the observed data... the points you would see on a plot
  int<lower=1> N1;
  vector[N1] x1;
  vector[N1] y1;

  # we predict values along x2. these would be the x values in between points.
  int<lower=1> N2;
  vector[N2] x2;
}
transformed data {
  int<lower=1> N;
  vector[N1 + N2] x;
  vector[N1 + N2] mu;
  N = N1 + N2;
  for (n in 1:N1) x[n] = x1[n];
  for (n in 1:N2) x[N1 + n] = x2[n];
  # fixed mean to 0
  for (i in 1:N) mu[i] = 0;
}
parameters {
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
  vector[N2] y2;
}
transformed parameters {
  real<lower=0> rho_sq;
  rho_sq = inv(inv_rho_sq);
}
model {
  matrix[N, N] Sigma;
  vector[N] y;

  eta_sq ~ normal(0, 5);
  inv_rho_sq ~ normal(0, 1);
  sigma_sq ~ normal(0, 5);

  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j], 2));
      Sigma[j, i] = Sigma[i, j];
    }
  }

  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = eta_sq + sigma_sq;

  for (n in 1:N1) y[n] = y1[n];
  for (n in 1:N2) y[N1 + n] = y2[n];

  y ~ multi_normal(mu, Sigma);
}

