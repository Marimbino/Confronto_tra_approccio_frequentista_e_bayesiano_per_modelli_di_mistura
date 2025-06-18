data {
  int<lower=1> N;    // Osservazioni
  int<lower=1> K;    // Componenti
  vector[N] y;       // Vettore osservazioni
}
transformed data {
  real y_mean = mean(y);   // media per prior
}
parameters {
  ordered[K] mu;   
  vector<lower=0>[K] sigma;
  simplex[K] p;   
}
model {
  mu ~ normal(y_mean, 3);   // Distribuzione normale 
  sigma ~ lognormal(1.3, 1);          // Distribuzione lognormale 
  p ~ dirichlet(rep_vector(1, K));    // Distribuzione di Dirichlet 
  for (i in 1:N) {
    vector[K] logf;
    for (j in 1:K)
      logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
    target += log_sum_exp(logf);
  }
}
