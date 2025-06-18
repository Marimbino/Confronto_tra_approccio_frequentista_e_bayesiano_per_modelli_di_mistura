data {
  int<lower=1> N;                  // Numero osservazioni
  int<lower=1> K;                  // Numero componenti
  matrix[N, 2] y;                  // Dati (bivariati)
}

transformed data {
  vector[2] y_mean;
  y_mean[1] = mean(col(y, 1));
  y_mean[2] = mean(col(y, 2));
}

parameters {
  array[K] vector[2] mu;                          // Medie
  array[K] vector<lower=0>[2] sigma;              // Deviazioni standard
  array[K] cholesky_factor_corr[2] L;             // Fattori di correlazione (Cholesky)
  simplex[K] p;                                   // Pesi della mistura
}

transformed parameters {
  array[K] matrix[2,2] L_Sigma;                   // Cholesky della matrice di covarianza
  for (k in 1:K)
    L_Sigma[k] = diag_pre_multiply(sigma[k], L[k]);
}

model {
  for (k in 1:K) {
    mu[k] ~ multi_normal(y_mean, diag_matrix(rep_vector(5.0, 2)));
    sigma[k] ~ lognormal(1, 0.7);
    L[k] ~ lkj_corr_cholesky(2);
  }
  p ~ dirichlet(rep_vector(2.0, K));

  for (n in 1:N) {
    vector[K] logf;
    for (k in 1:K)
      logf[k] = log(p[k]) + multi_normal_cholesky_lpdf(y[n] | mu[k], L_Sigma[k]);
    target += log_sum_exp(logf);
  }
}

generated quantities {
  array[K] cov_matrix[2] Sigma;
  for (k in 1:K)
    Sigma[k] = L_Sigma[k] * L_Sigma[k]';
}
