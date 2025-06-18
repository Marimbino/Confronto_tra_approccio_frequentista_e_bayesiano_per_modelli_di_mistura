data {
  int<lower=1> N;           
  int<lower=1> K;            
  matrix[N, 2] y;             
}
transformed data {
  vector[2] y_mean;
  y_mean[1] = mean(col(y, 1));
  y_mean[2] = mean(col(y, 2));
}
parameters {
  array[K] vector[2] mu;                
  array[K] cov_matrix[2] sigma;         
  simplex[K] p;                         // Pesi mistura
}
model {
  for (k in 1:K) {
    mu[k] ~ multi_normal(y_mean, diag_matrix(rep_vector(5, 2)));    // Prior sulla media
    sigma[k] ~ inv_wishart(10, diag_matrix(rep_vector(3, 2)));       // Prior debole su covarianza
  }

  p ~ dirichlet(rep_vector(2, K));    // Prior moderatamente informativa

  for (n in 1:N) {
    vector[K] logf;
    for (k in 1:K) {
      logf[k] = log(p[k]) + multi_normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(logf);
  }
}
