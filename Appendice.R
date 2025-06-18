# Confronto tra approccio frequentista e bayesiano per modelli di mistura

library(mixtools)
library(mclust)
library(cmdstanr)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

R <- 100

# Funzione per simulare R volte
sim_2comp <- function(N, mu.true, sigma.true, p.true, cmd.mod) {
   
   # Tibble risultati
   ris <- tibble(R = integer(), Metodo = factor(),
                  mu1 = numeric(), mu2 = numeric(),
                  sigma1 = numeric(), sigma2 = numeric(),
                  p1 = numeric(), p2 = numeric(),
                  tempo = numeric())
   
   set.seed(123)
   for (i in 1:R) {
      
      # Labels e valori
      z <- rbinom(N, 1, p.true[2]) + 1 
      y <- rnorm(N, mean = mu.true[z], sd = sigma.true[z]) 
      
      # Approccio classico - Algoritmo EM (normalmixEM)
      start <- Sys.time()
      fit.em <- normalmixEM(x = y, k = 2, maxit = 1000, epsilon = 1e-8)
      tempo.em <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      ris <- add_row(ris, R = i, Metodo = "EM",
                     mu1 = min(fit.em$mu),
                     mu2 = max(fit.em$mu),
                     sigma1 = fit.em$sigma[which.min(fit.em$mu)],
                     sigma2 = fit.em$sigma[which.max(fit.em$mu)],
                     p1 = fit.em$lambda[which.min(fit.em$mu)],
                     p2 = fit.em$lambda[which.max(fit.em$mu)],
                     tempo = tempo.em)
      
      # Approccio bayesiano - Metodo MCMC (cmdstan_model$sample)
      stan.data <- list(N = N, K = 2, y = y)
      start <- Sys.time()
      fit.mcmc <- cmd.mod$sample(data = stan.data,
                                 iter_sampling = 2000, iter_warmup = 1000,
                                 chains = 4, parallel_chains = 4,
                                 seed = 123)
      tempo.mcmc <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.mcmc.sum <- fit.mcmc$summary()
      ris <- add_row(ris, R = i, Metodo = "MCMC",
                     mu1 = as.numeric(fit.mcmc.sum[2, 2]),
                     mu2 = as.numeric(fit.mcmc.sum[3, 2]),
                     sigma1 = as.numeric(fit.mcmc.sum[4, 2]),
                     sigma2 = as.numeric(fit.mcmc.sum[5, 2]),
                     p1 = as.numeric(fit.mcmc.sum[6, 2]),
                     p2 = as.numeric(fit.mcmc.sum[7, 2]),
                     tempo = tempo.mcmc)
   }
   return(ris)
}

# Funzione per ottenere le metriche
met_2comp <- function(ris, mu.true, sigma.true, p.true){
   
   ris.med <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")), 
                   names_to = "Parametro", values_to = "Stima") |>
      group_by(Metodo, Parametro) |>
      summarize(Media = mean(Stima)) |>
      arrange(Metodo, Parametro)
   
   ris.mse <- ris |>
      mutate(bias.mu1 = mu1 - mu.true[1], bias.mu2 = mu2 - mu.true[2],
             bias.sigma1 = sigma1 - sigma.true[1], 
             bias.sigma2 = sigma2 - sigma.true[2],
             bias.p1 = p1 - p.true[1], bias.p2 = p2 - p.true[2]) |>
      pivot_longer(cols = starts_with("bias."), names_to = "Parametro", 
                   names_prefix = "bias.", values_to = "Errore") |>
      group_by(Parametro, Metodo) |>
      summarize(BIAS = mean(Errore), SD = sd(Errore),
                MSE = mean(Errore^2), .groups = "drop")
   
   ris.tot <- ris.med |> 
      full_join(ris.mse)  |>
      arrange(Parametro, Metodo)
   
   return(ris.tot)
}

# Funzione per boxplot
plot_2comp <- function(ris, mu.true, sigma.true, p.true){
   
   par.true <- tibble(Parametro = c("mu1", "mu2", "sigma1", "sigma2", "p1", "p2"),
                      val = c(mu.true, sigma.true, p.true))
   
   ris.long <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")),
                   names_to = "Parametro", values_to = "Stima")
   
   gg <- ggplot(ris.long, aes(x = Metodo, y = Stima, fill = Metodo)) +
      geom_boxplot(width = 0.3, alpha = 0.6, outliers = F) +
      geom_hline(data = par.true, aes(yintercept = val),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
      facet_wrap(~ Parametro, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("EM" = "#e41a1c", "MCMC" = "#377eb8")) +
      labs(title = "Confronto stime EM e MCMC", y = "Stima", x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
   
   return(gg)
}


# Simulazione 1 -----------------------------------------------------------

# Codice stan (sim1.stan)
data {
   int<lower=1> N;
   int<lower=1> K;
   vector[N] y;
}
transformed data {
   real y_mean = mean(y);  
}
parameters {
   ordered[K] mu;   
   vector<lower=0>[K] sigma;
   simplex[K] p;   
}
model {
   mu ~ normal(y_mean, 2);   
   sigma ~ lognormal(0, 0.5);       
   p ~ dirichlet(rep_vector(2, K));    
   for (i in 1:N) {
      vector[K] logf;
      for (j in 1:K)
         logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 500
p.true <- c(0.5, 0.5) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 1)
cmdstan.mod1 <- cmdstan_model("sim1.stan")

ris1 <- read.csv("sim1_ris.csv", header = T, sep = ",")
ris1 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris1, mu.true, sigma.true, p.true)
ris1 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris1, mu.true, sigma.true, p.true)
rm(ris1)


# Simulazione 2 -----------------------------------------------------------

N <- 100
p.true <- c(0.5, 0.5) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 1)

ris2 <- read.csv("sim2_ris.csv", header = T, sep = ",")
ris2 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris2, mu.true, sigma.true, p.true)
ris2 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris2, mu.true, sigma.true, p.true)
rm(ris2)


# Simulazione 3 -----------------------------------------------------------

N <- 500
p.true <- c(0.5, 0.5) 
mu.true <- c(4, 6) 
sigma.true <- c(1, 1)

ris3 <- read.csv("sim3_ris.csv", header = T, sep = ",")
ris3 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris3, mu.true, sigma.true, p.true)
ris3 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris3, mu.true, sigma.true, p.true)
rm(ris3)


# Simulazione 4 -----------------------------------------------------------

N <- 100
p.true <- c(0.5, 0.5) 
mu.true <- c(4, 6) 
sigma.true <- c(1, 1)

ris4 <- read.csv("sim4_ris.csv", header = T, sep = ",")
ris4 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris4, mu.true, sigma.true, p.true)
ris4 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris4, mu.true, sigma.true, p.true)
rm(ris4)


# Simulazione 5 -----------------------------------------------------------

N <- 500
p.true <- c(0.85, 0.15) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 1)

ris5 <- read.csv("sim5_ris.csv", header = T, sep = ",")
ris5 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris5, mu.true, sigma.true, p.true)
ris5 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris5, mu.true, sigma.true, p.true)
rm(ris5)


# Simulazione 6 -----------------------------------------------------------

# Codice stan (sim6.stan)
data {
   int<lower=1> N;
   int<lower=1> K;
   vector[N] y;
}
transformed data {
   real y_mean = mean(y);
}
parameters {
   ordered[K] mu;   
   vector<lower=0>[K] sigma;
   simplex[K] p;   
}
model {
   mu ~ normal(y_mean, 2);
   sigma ~ lognormal(0, 0.5);
   p ~ dirichlet(rep_vector(1, K));
   for (i in 1:N) {
      vector[K] logf;
      for (j in 1:K)
         logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 100
p.true <- c(0.85, 0.15) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 1)
cmdstan.mod6 <- cmdstan_model("sim6.stan")

ris6 <- read.csv("sim6_ris.csv", header = T, sep = ",")
ris6 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod6)
met_2comp(ris6, mu.true, sigma.true, p.true)
ris6 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris6, mu.true, sigma.true, p.true)
rm(ris6)


# Simulazione 7 -----------------------------------------------------------

N <- 500
p.true <- c(0.5, 0.5) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 3)

ris7 <- read.csv("sim7_ris.csv", header = T, sep = ",")
ris7 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod1)
met_2comp(ris7, mu.true, sigma.true, p.true)
ris7 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris7, mu.true, sigma.true, p.true)
rm(ris7)


# Simulazione 8 -----------------------------------------------------------

# Codice stan (sim8.stan)
data {
   int<lower=1> N;
   int<lower=1> K;
   vector[N] y;
}
transformed data {
   real y_mean = mean(y); 
}
parameters {
   ordered[K] mu;   
   vector<lower=0>[K] sigma;
   simplex[K] p;   
}
model {
   mu ~ normal(y_mean, 3);
   sigma ~ lognormal(1.3, 1);
   p ~ dirichlet(rep_vector(1, K));
   for (i in 1:N) {
      vector[K] logf;
      for (j in 1:K)
         logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 500
p.true <- c(0.5, 0.5) 
mu.true <- c(2, 6) 
sigma.true <- c(1, 1)
cmdstan.mod8 <- cmdstan_model("sim8.stan")

ris8 <- read.csv("sim8_ris.csv", header = T, sep = ",")
ris8 <- sim_2comp(N, mu.true, sigma.true, p.true, cmdstan.mod8)
met_2comp(ris8, mu.true, sigma.true, p.true)
ris8 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_2comp(ris8, mu.true, sigma.true, p.true)
rm(ris8)


# Simulazione 9 -----------------------------------------------------------

# Funzione per simulare R volte
sim_3comp <- function(N, mu.true, sigma.true, p.true, cmd.mod) {
   
   # Tibble risultati
   ris <- tibble(R = integer(), Metodo = factor(),
                 mu1 = numeric(), mu2 = numeric(), mu3 = numeric(),
                 sigma1 = numeric(), sigma2 = numeric(), sigma3 = numeric(),
                 p1 = numeric(), p2 = numeric(), p3 = numeric(),
                 tempo = numeric())
   
   set.seed(123)
   for (i in 1:R) {
      
      # Labels e valori
      z <- sample(1:3, size=N, prob = p.true, replace=T)
      y <- rnorm(N, mean = mu.true[z], sd = sigma.true[z]) 
      
      # Approccio classico - Algoritmo EM (normalmixEM)
      start <- Sys.time()
      fit.em <- normalmixEM(x = y, k = 3, maxit = 1000, epsilon = 1e-8)
      tempo.em <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.em.sort <- matrix(c(fit.em$mu, fit.em$sigma, fit.em$lambda), 
                            ncol = 3, byrow = T)
      fit.em.sort <- fit.em.sort[, order(fit.em.sort[1, ])]
      ris <- add_row(ris, R = i, Metodo = "EM",
                     mu1 = fit.em.sort[1, 1], mu2 = fit.em.sort[1, 2],
                     mu3 = fit.em.sort[1, 3], sigma1 = fit.em.sort[2, 1], 
                     sigma2 = fit.em.sort[2, 2], sigma3 = fit.em.sort[2, 3],
                     p1 = fit.em.sort[3, 1], p2 = fit.em.sort[3, 2],
                     p3 = fit.em.sort[3, 3], tempo = tempo.em)
      
      # Approccio bayesiano - Metodo MCMC (cmdstan_model$sample)
      stan.data <- list(N = N, K = 3, y = y)
      start <- Sys.time()
      fit.mcmc <- cmd.mod$sample(data = stan.data,
                                 iter_sampling = 2000, iter_warmup = 1000,
                                 chains = 4, parallel_chains = 4,
                                 seed = 123)
      tempo.mcmc <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.mcmc.sort <- matrix(fit.mcmc$summary()$mean[-1], ncol = 3, byrow = T)
      fit.mcmc.sort <- fit.mcmc.sort[, order(fit.mcmc.sort[1, ])]
      ris <- add_row(ris, R = i, Metodo = "MCMC",
                      mu1 = fit.mcmc.sort[1, 1], mu2 = fit.mcmc.sort[1, 2],
                      mu3 = fit.mcmc.sort[1, 3], sigma1 = fit.mcmc.sort[2, 1],
                      sigma2 = fit.mcmc.sort[2, 2], sigma3 = fit.mcmc.sort[2, 3],
                      p1 = fit.mcmc.sort[3, 1], p2 = fit.mcmc.sort[3, 2],
                      p3 = fit.mcmc.sort[3, 3], tempo = tempo.mcmc)
   }
   return(ris)
}

# Funzione per ottenere le metriche
met_3comp <- function(ris, mu.true, sigma.true, p.true){
   
   ris.med <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")), 
                   names_to = "Parametro", values_to = "Stima") |>
      group_by(Metodo, Parametro) |>
      summarize(Media = mean(Stima)) |>
      arrange(Metodo, Parametro)
   
   ris.mse <- ris |>
      mutate(bias.mu1 = mu1 - mu.true[1], bias.mu2 = mu2 - mu.true[2],
             bias.mu3 = mu3 - mu.true[3], bias.sigma1 = sigma1 - sigma.true[1], 
             bias.sigma2 = sigma2 - sigma.true[2],
             bias.sigma3 = sigma3 - sigma.true[3],
             bias.p1 = p1 - p.true[1], bias.p2 = p2 - p.true[2],
             bias.p3 = p3 - p.true[3]) |>
      pivot_longer(cols = starts_with("bias."), names_to = "Parametro", 
                   names_prefix = "bias.", values_to = "Errore") |>
      group_by(Parametro, Metodo) |>
      summarize(BIAS = mean(Errore), SD = sd(Errore),
                MSE = mean(Errore^2), .groups = "drop")
   
   ris.tot <- ris.med |> 
      full_join(ris.mse)  |>
      arrange(Parametro, Metodo)
   
   return(ris.tot)
}

# Funzione per boxplot
plot_3comp <- function(ris, mu.true, sigma.true, p.true){
   
   par.true <- tibble(Parametro = c("mu1", "mu2", "mu3", "sigma1",
                                    "sigma2", "sigma3", "p1", "p2", "p3"),
                      val = c(mu.true, sigma.true, p.true))
   
   ris.long <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")),
                   names_to = "Parametro", values_to = "Stima")
   
   gg <- ggplot(ris.long, aes(x = Metodo, y = Stima, fill = Metodo)) +
      geom_boxplot(width = 0.3, alpha = 0.6, outliers = F) +
      geom_hline(data = par.true, aes(yintercept = val),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
      facet_wrap(~ Parametro, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("EM" = "#e41a1c", "MCMC" = "#377eb8")) +
      labs(title = "Confronto stime EM e MCMC", y = "Stima", x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
   
   return(gg)
}

# Codice stan (sim9.stan)
data {
   int<lower=1> N;
   int<lower=1> K;
   vector[N] y;
}
transformed data {
   real y_mean = mean(y);
}
parameters {
   ordered[K] mu;   
   vector<lower=0>[K] sigma;
   simplex[K] p;   
}
model {
   mu ~ normal(y_mean, 3); 
   sigma ~ lognormal(0, 0.5); 
   p ~ dirichlet(rep_vector(2, K)); 
   for (i in 1:N) {
      vector[K] logf;
      for (j in 1:K)
         logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 750
p.true <- c(1, 1, 1)/3 
mu.true <- c(1, 4, 7)
sigma.true <- c(1, 1, 1)
cmdstan.mod9 <- cmdstan_model("sim9.stan")

ris9 <- read.csv("sim9_ris.csv", header = T, sep = ",")
ris9 <- sim_3comp(N, mu.true, sigma.true, p.true, cmdstan.mod9)
met_3comp(ris9, mu.true, sigma.true, p.true)
ris9 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_3comp(ris9, mu.true, sigma.true, p.true)
rm(ris9)


# Simulazione 10 ----------------------------------------------------------

# Funzione per simulare R volte
sim_3comp_a <- function(N, mu.true, sigma.true, p.true, cmd.mod) {
   
   # Tibble risultati
   ris <- tibble(R = integer(), Metodo = factor(),
                 mu1 = numeric(), mu2 = numeric(), mu3 = numeric(),
                 sigma1 = numeric(), sigma2 = numeric(), sigma3 = numeric(),
                 p1 = numeric(), p2 = numeric(), p3 = numeric(),
                 tempo = numeric())
   
   set.seed(123)
   for (i in 1:R) {
      
      # Labels e valori
      z <- sample(1:3, size=N, prob = p.true, replace=T)
      y <- rnorm(N, mean = mu.true[z], sd = sigma.true[z]) 
      
      km <- kmeans(x = y, centers = 3)
      inizial <- matrix(c(km$centers, km$size / N), ncol = 3, byrow = T)
      inizial <- inizial[, order(inizial[1, ])]
      
      # Approccio classico - Algoritmo EM (normalmixEM)
      start <- Sys.time()
      fit.em <- normalmixEM(x = y, k = 3, maxit = 1000, epsilon = 1e-8)
      tempo.em <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.em.sort <- matrix(c(fit.em$mu, fit.em$sigma, fit.em$lambda), 
                            ncol = 3, byrow = T)
      fit.em.sort <- fit.em.sort[, order(fit.em.sort[1, ])]
      ris <- add_row(ris, R = i, Metodo = "EM",
                     mu1 = fit.em.sort[1, 1], mu2 = fit.em.sort[1, 2],
                     mu3 = fit.em.sort[1, 3], sigma1 = fit.em.sort[2, 1], 
                     sigma2 = fit.em.sort[2, 2], sigma3 = fit.em.sort[2, 3],
                     p1 = fit.em.sort[3, 1], p2 = fit.em.sort[3, 2],
                     p3 = fit.em.sort[3, 3], tempo = tempo.em)
      
      # Approccio bayesiano - Metodo MCMC (cmdstan_model$sample)
      stan.data <- list(N = N, K = 3, y = y)
      start <- Sys.time()
      fit.mcmc <- cmd.mod$sample(data = stan.data,
                                 iter_sampling = 2000, iter_warmup = 1000,
                                 chains = 4, parallel_chains = 4,
                                 seed = 123)
      tempo.mcmc <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.mcmc.sort <- matrix(fit.mcmc$summary()$mean[-1], ncol = 3, byrow = T)
      fit.mcmc.sort <- fit.mcmc.sort[, order(fit.mcmc.sort[1, ])]
      ris <- add_row(ris, R = i, Metodo = "MCMC",
                     mu1 = fit.mcmc.sort[1, 1], mu2 = fit.mcmc.sort[1, 2],
                     mu3 = fit.mcmc.sort[1, 3], sigma1 = fit.mcmc.sort[2, 1],
                     sigma2 = fit.mcmc.sort[2, 2], sigma3 = fit.mcmc.sort[2, 3],
                     p1 = fit.mcmc.sort[3, 1], p2 = fit.mcmc.sort[3, 2],
                     p3 = fit.mcmc.sort[3, 3], tempo = tempo.mcmc)
   }
   return(ris)
}

# Codice stan (sim10.stan)
data {
   int<lower=1> N;
   int<lower=1> K;
   vector[N] y;
   vector[K] mu_km;
   vector[K] p_km;
}
transformed data {
   real y_mean = mean(y);
}
parameters {
   ordered[K] mu;   
   vector<lower=0>[K] sigma;
   simplex[K] p;   
}
model {
   mu ~ normal(mu_km, 1.5); 
   sigma ~ lognormal(1, 0.8); 
   p ~ dirichlet(p_km); 
   for (i in 1:N) {
      vector[K] logf;
      for (j in 1:K)
         logf[j] = log(p[j]) + normal_lpdf(y[i] | mu[j], sigma[j]);
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 300
p.true <- c(0.15, 0.5, 0.35)
mu.true <- c(1, 6, 7)
sigma.true <- c(3, 1, 0.5)
cmdstan.mod10 <- cmdstan_model("sim10.stan")

ris10 <- read.csv("sim10_ris.csv", header = T, sep = ",")
ris10 <- sim_3comp(N, mu.true, sigma.true, p.true, cmdstan.mod10)
met_3comp(ris10, mu.true, sigma.true, p.true)
ris10 |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_3comp(ris10, mu.true, sigma.true, p.true)
rm(ris10)


# Simulazione 11 ----------------------------------------------------------

# Funzione per campionare da normale bivariata
r_mixNorm_mvt <- function(n, mu1, mu2, mu3, sigma1, sigma2, sigma3, p){
   
   componenti <- sample(1:3, size=n, prob = p, replace=T)
   n1 <- length(which(componenti==1)) #n * p1
   n2 <- length(which(componenti==2)) #n * p2
   n3 <- length(which(componenti==3)) #n * (1-p1-p2)
   
   y <- matrix(NA, nrow=n, ncol=2)
   y[which(componenti==1),1:2] <- rmvnorm(n1, mu1, sigma1)   
   y[which(componenti==2),1:2] <- rmvnorm(n2, mu2, sigma2)
   y[which(componenti==3),1:2] <- rmvnorm(n3, mu3, sigma3)
   return(y)
}

# Funzione per simulare R volte
sim_3comp_2d <- function(N, mu.true, sigma.true, p.true, cmd.mod) {
   
   # Tibble risultati
   ris11 <- tibble(R = integer(), Metodo = factor(),
                   mu1a = numeric(), mu1b = numeric(), 
                   mu2a = numeric(), mu2b = numeric(),
                   mu3a = numeric(), mu3b = numeric(),
                   sigma1a = numeric(), sigma1b = numeric(), 
                   sigma1c = numeric(), sigma1d = numeric(),
                   sigma2a = numeric(), sigma2b = numeric(), 
                   sigma2c = numeric(), sigma2d = numeric(),
                   sigma3a = numeric(), sigma3b = numeric(), 
                   sigma3c = numeric(), sigma3d = numeric(),
                   p1 = numeric(), p2 = numeric(), p3 = numeric(),
                   tempo = numeric()
   )
   
   set.seed(123)
   for (i in 1:R) {
      
      # Valori
      y <- r_mixNorm_mvt(N, mu1.true, mu2.true, mu3.true, 
                         sigma1.true, sigma2.true, sigma3.true, p.true) 
      
      # Approccio classico - Algoritmo EM (normalmixEM)
      start <- Sys.time()
      fit.em <- fit.em <- Mclust(data = y, G = 3)
      tempo.em <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.em <- fit.em$parameters
      
      fit.em.m <- matrix(fit.em$mean, ncol = 2, byrow = T)
      fit.em.v <- matrix(fit.em$variance$sigma, ncol = 4, byrow = T)
      fit.em.p <- matrix(fit.em$pro, ncol = 1)
      
      fit.em.sort <- matrix(cbind(fit.em.m, fit.em.v, fit.em.p), nrow = 3)
      fit.em.sort <- fit.em.sort[order(fit.em.sort[,1]), ]
      
      ris <- add_row(ris, R = i, Metodo = "EM",
                     mu1a = fit.em.sort[1, 1], mu1b = fit.em.sort[1, 2],
                     mu2a = fit.em.sort[2, 1], mu2b = fit.em.sort[2, 2],
                     mu3a = fit.em.sort[3, 1], mu3b = fit.em.sort[3, 2],
                     sigma1a = fit.em.sort[1, 3], sigma1b = fit.em.sort[1, 4],
                     sigma1c = fit.em.sort[1, 5], sigma1d = fit.em.sort[1, 6],
                     sigma2a = fit.em.sort[2, 3], sigma2b = fit.em.sort[2, 4],
                     sigma2c = fit.em.sort[2, 5], sigma2d = fit.em.sort[2, 6],
                     sigma3a = fit.em.sort[3, 3], sigma3b = fit.em.sort[3, 4],
                     sigma3c = fit.em.sort[3, 5], sigma3d = fit.em.sort[3, 6],
                     p1 = fit.em.sort[1, 7], p2 = fit.em.sort[2, 7],
                     p3 = fit.em.sort[3, 7], tempo = tempo.em)
      
      # Approccio bayesiano - Metodo MCMC (cmdstan_model$sample)
      stan.data <- list(N = N, K = 3, y = y)
      start <- Sys.time()
      fit.mcmc <- cmd.mod$sample(data = stan.data,
                                 iter_sampling = 2000, iter_warmup = 1000,
                                 chains = 4, parallel_chains = 4,
                                 seed = 123)
      tempo.mcmc <- as.numeric(difftime(Sys.time(), start, units="secs"))
      
      fit.mcmc <- fit.mcmc$summary()$mean[-1]
      
      fit.mcmc.m <- matrix(fit.mcmc[1:6], ncol = 2, byrow = F)
      fit.mcmc.v <- matrix(fit.mcmc[7:18], ncol = 4, byrow = F)
      fit.mcmc.p <- matrix(fit.mcmc[19:21], ncol = 1)
      
      fit.mcmc.sort <- matrix(cbind(fit.mcmc.m, fit.mcmc.v, fit.mcmc.p), nrow = 3)
      fit.mcmc.sort <- fit.mcmc.sort[order(fit.mcmc.sort[,1]), ]
      
      ris <- add_row(ris, R = i, Metodo = "MCMC",
                     mu1a = fit.mcmc.sort[1, 1], mu1b = fit.mcmc.sort[1, 2],
                     mu2a = fit.mcmc.sort[2, 1], mu2b = fit.mcmc.sort[2, 2],
                     mu3a = fit.mcmc.sort[3, 1], mu3b = fit.mcmc.sort[3, 2],
                     sigma1a = fit.mcmc.sort[1, 3], sigma1b = fit.mcmc.sort[1, 4],
                     sigma1c = fit.mcmc.sort[1, 5], sigma1d = fit.mcmc.sort[1, 6],
                     sigma2a = fit.mcmc.sort[2, 3], sigma2b = fit.mcmc.sort[2, 4],
                     sigma2c = fit.mcmc.sort[2, 5], sigma2d = fit.mcmc.sort[2, 6],
                     sigma3a = fit.mcmc.sort[3, 3], sigma3b = fit.mcmc.sort[3, 4],
                     sigma3c = fit.mcmc.sort[3, 5], sigma3d = fit.mcmc.sort[3, 6],
                     p1 = fit.mcmc.sort[1, 7], p2 = fit.mcmc.sort[2, 7],
                     p3 = fit.mcmc.sort[3, 7],tempo = tempo.mcmc)
   }
   return(ris)
}

# Funzione per ottenere le metriche
met_3comp_2d <- function(ris, mu.true, sigma.true, p.true){
   
   ris.med <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")), 
                   names_to = "Parametro", values_to = "Stima") |>
      group_by(Metodo, Parametro) |>
      summarize(Media = mean(Stima)) |>
      arrange(Metodo, Parametro)
   
   ris.mse <- ris |>
      mutate(bias.mu1a = mu1a - mu1.true[1], bias.mu1b = mu1b - mu1.true[2],
             bias.mu2a = mu2a - mu2.true[1], bias.mu2b = mu2b - mu2.true[2],
             bias.mu3a = mu3a - mu3.true[1], bias.mu3b = mu3b - mu3.true[2],
             bias.sigma1a = sigma1a - sigma1.true[1], 
             bias.sigma1b = sigma1b - sigma1.true[2], 
             bias.sigma1c = sigma1c - sigma1.true[3], 
             bias.sigma1d = sigma1d - sigma1.true[4],
             bias.sigma2a = sigma2a - sigma2.true[1], 
             bias.sigma2b = sigma2b - sigma2.true[2], 
             bias.sigma2c = sigma2c - sigma2.true[3], 
             bias.sigma2d = sigma2d - sigma2.true[4],
             bias.sigma3a = sigma3a - sigma3.true[1], 
             bias.sigma3b = sigma3b - sigma3.true[2], 
             bias.sigma3c = sigma3c - sigma3.true[3], 
             bias.sigma3d = sigma3d - sigma3.true[4],
             bias.p1 = p1 - p.true[1], bias.p2 = p2 - p.true[2],
             bias.p3 = p3 - p.true[3]) |>
      pivot_longer(cols = starts_with("bias."), names_to = "Parametro", 
                   names_prefix = "bias.", values_to = "Errore") |>
      group_by(Parametro, Metodo) |>
      summarize(BIAS = mean(Errore), SD = sd(Errore),
                MSE = mean(Errore^2), .groups = "drop")
   
   ris.tot <- ris.med |> 
      full_join(ris.mse)  |>
      arrange(Parametro, Metodo)
   
   return(ris.tot)
}

# Funzione per boxplot
plot_3comp_2d <- function(ris, mu.true, sigma.true, p.true){
   
   par.true <- tibble(Parametro = c("mu1a", "mu1b", "mu2a", "mu2b", "mu3a", "mu3b",
                                    "sigma1a", "sigma1b", "sigma1c", "sigma1d",
                                    "sigma2a", "sigma2b", "sigma2c", "sigma2d",
                                    "sigma3a", "sigma3b", "sigma3c", "sigma3d",
                                    "p1", "p2", "p3"),
                      val = c(mu1.true, mu2.true, mu3.true,
                              sigma1.true, sigma2.true, sigma3.true, p.true))
   
   ris.long <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")),
                   names_to = "Parametro", values_to = "Stima")
   
   gg <- ggplot(ris.long, aes(x = Metodo, y = Stima, fill = Metodo)) +
      geom_boxplot(width = 0.3, alpha = 0.6, outliers = F) +
      geom_hline(data = par.true, aes(yintercept = val),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
      facet_wrap(~ Parametro, scales = "free_y", ncol = 7) +
      scale_fill_manual(values = c("EM" = "#e41a1c", "MCMC" = "#377eb8")) +
      labs(title = "Confronto stime EM e MCMC", y = "Stima", x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
   
   return(gg)
}

# Codice stan (sim11.stan)
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
   simplex[K] p;
}
model {
   for (k in 1:K) {
      mu[k] ~ multi_normal(y_mean, diag_matrix(rep_vector(5, 2)));
      sigma[k] ~ inv_wishart(10, diag_matrix(rep_vector(3, 2)));
   }
   
   p ~ dirichlet(rep_vector(2, K));
   
   for (n in 1:N) {
      vector[K] logf;
      for (k in 1:K) {
         logf[k] = log(p[k]) + multi_normal_lpdf(y[n] | mu[k], sigma[k]);
      }
      target += log_sum_exp(logf);
   }
}

# Codice R
N <- 300
p.true <- rep(1/3, 3) 
mu1.true <- c(1, 5)
mu2.true <- c(4, 1)
mu3.true <- c(6, 7)
sigma1.true <- matrix(c(2,2,2,6), ncol=2)
sigma2.true <- matrix(c(3,1,1,1), ncol=2)
sigma3.true <- matrix(c(3,-3,-3,5), ncol=2)
cmdstan.mod11 <- cmdstan_model("sim11.stan")

ris11 <- read.csv("sim11_ris.csv", header = T, sep = ",")
ris11 <- sim_3comp(N, mu.true, sigma.true, p.true, cmdstan.mod11)
print(met_3comp_2d(ris11, mu.true, sigma.true, p.true), n=50)
ris11 |> 
   filter(tempo < 1000) |> 
   group_by(Metodo) |> 
   summarize(mean(tempo))
plot_3comp_2d(ris11, mu.true, sigma.true, p.true)
rm(ris11)


# Simulazione 12 ----------------------------------------------------------

# Codice stan (sim12.stan)
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
   array[K] vector<lower=0>[2] sigma;
   array[K] cholesky_factor_corr[2] L;
   simplex[K] p;
}
transformed parameters {
   array[K] matrix[2,2] L_Sigma;
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


# Codice R'
   N <- 300
   p.true <- rep(1/3, 3) 
   mu1.true <- c(1, 5)
   mu2.true <- c(4, 1)
   mu3.true <- c(6, 7)
   sigma1.true <- matrix(c(2,2,2,6), ncol=2)
   sigma2.true <- matrix(c(3,1,1,1), ncol=2)
   sigma3.true <- matrix(c(3,-3,-3,5), ncol=2)
   cmdstan.mod12 <- cmdstan_model("sim12.stan")
   
   ris12 <- read.csv("sim12_ris.csv", header = T, sep = ",")
   ris12 <- sim_3comp(N, mu.true, sigma.true, p.true, cmdstan.mod12)
   print(met_3comp_2d(ris12, mu.true, sigma.true, p.true), n=50)
   ris12 |> 
      filter(tempo < 1000) |>
      group_by(Metodo) |> 
      summarize(mean(tempo))
   plot_3comp_2d(ris12, mu.true, sigma.true, p.true)
   rm(ris12)
   

# Grafico boxplot sim11 e sim12
plot_3comp_2d.3 <- function(ris.a, ris.b, mu.true, sigma.true, p.true){
   
   par.true <- tibble(Parametro = c("mu1a", "mu1b", "mu2a", "mu2b", "mu3a", "mu3b",
                                    "sigma1a", "sigma1b", "sigma1c", "sigma1d",
                                    "sigma2a", "sigma2b", "sigma2c", "sigma2d",
                                    "sigma3a", "sigma3b", "sigma3c", "sigma3d",
                                    "p1", "p2", "p3"),
                      val = c(mu1.true, mu2.true, mu3.true,
                              sigma1.true, sigma2.true, sigma3.true, p.true))
   
   ris.a <- ris11 |> 
      mutate(Metodo = if_else(Metodo == "MCMC", "WIS", Metodo))
   ris.b <- ris12 |> 
      filter(Metodo == "MCMC") |> 
      mutate(Metodo = "LKJ")
   ris <- ris.a |> 
      add_row(ris.b)
   
   ris.long <- ris |>
      pivot_longer(cols = starts_with(c("mu", "sigma", "p")),
                   names_to = "Parametro", values_to = "Stima")
   
   gg <- ggplot(ris.long, aes(x = Metodo, y = Stima, fill = Metodo)) +
      geom_boxplot(width = 0.3, alpha = 0.6, outliers = F) +
      geom_hline(data = par.true, aes(yintercept = val),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
      facet_wrap(~ Parametro, scales = "free_y", ncol = 7) +
      scale_fill_manual(values = c("EM" = "#e41a1c", "WIS" = "#377eb8", 
                                   "LKJ" = "#4daf4a")) +
      labs(title = "Confronto stime EM, Invert-Wishart e LKJ", y = "Stima", x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
   
   return(gg)
}

plot_3comp_2d.3(ris11, ris12, mu.true, sigma.true, p.true)