# Funzioni base
r_mixNorm_1 <- function(n, mu1, mu2, sd1, sd2, p){
   
   componenti <- sample(1:2, size=n, prob = p, replace=T)
   n1 <- length(which(componenti==1)) #n * p1
   n2 <- length(which(componenti==2)) #n * (1-p1)
   
   y <- matrix(NA, nrow=n, ncol=2)
   y[which(componenti==1),1] <- rnorm(n1, mu1, sd1)   
   y[which(componenti==2),1] <- rnorm(n2, mu2, sd2)
   y[,2] <- componenti   
   return(y)
}


# Simulazione 1 e 8 -------------------------------------------------------

n <- 500
p <- c(0.5, 0.5)
mu1 <- 2
mu2 <- 6
sd1 <- 1
sd2 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
   geom_histogram(aes(y = after_stat(density)), 
                  fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
   geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
   labs(title = "Simulazione 1", x = "y", y = "Densità") +
   theme_minimal(base_size = 13)
ggsave("sim1_densità.pdf", width = 7, height = 4)



# Simulazione 2 -----------------------------------------------------------

n <- 100
p <- c(0.5, 0.5)
mu1 <- 2
mu2 <- 6
sd1 <- 1
sd2 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
   geom_histogram(aes(y = after_stat(density)), 
                  fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
   geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
   labs(title = "Simulazione 2", x = "y", y = "Densità") +
   theme_minimal(base_size = 13)


# Simulazione 3 -----------------------------------------------------------

n <- 500
p <- c(0.5, 0.5)
mu1 <- 4
mu2 <- 6
sd1 <- 1
sd2 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 3", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)


# Simulazione 4 -----------------------------------------------------------

n <- 100
p <- c(0.5, 0.5)
mu1 <- 4
mu2 <- 6
sd1 <- 1
sd2 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 4", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)


# Simulazione 5 e 6 -------------------------------------------------------

n <- 500
p <- c(0.85, 0.15)
mu1 <- 2
mu2 <- 6
sd1 <- 1
sd2 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 5", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)


# Simulazione 7 -----------------------------------------------------------

n <- 500
p <- c(0.5, 0.5)
mu1 <- 2
mu2 <- 6
sd1 <- 1
sd2 <- 3

set.seed(123)
y <- as.data.frame(r_mixNorm_1(n, mu1, mu2, sd1, sd2, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 7", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)



# 3 componenti ------------------------------------------------------------

r_mixNorm_2 <- function(n, mu1, mu2, mu3, sd1, sd2, sd3, p){
   
   componenti <- sample(1:3, size=n, prob = p, replace=T)
   n1 <- length(which(componenti==1)) #n * p1
   n2 <- length(which(componenti==2)) #n * p2
   n3 <- length(which(componenti==3)) #n * p3
   
   y <- matrix(NA, nrow=n, ncol=2)
   y[which(componenti==1),1] <- rnorm(n1, mu1, sd1)   
   y[which(componenti==2),1] <- rnorm(n2, mu2, sd2)
   y[which(componenti==3),1] <- rnorm(n3, mu3, sd3)
   y[,2] <- componenti   
   return(y)
}


# Simulazione 9 -----------------------------------------------------------

n <- 750
p <- c(1, 1, 1)/3
mu1 <- 1
mu2 <- 4
mu3 <- 7
sd1 <- 1
sd2 <- 1
sd3 <- 1

set.seed(123)
y <- as.data.frame(r_mixNorm_2(n, mu1, mu2, mu3, sd1, sd2, sd3, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 9", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)


# Simulazione 10 ----------------------------------------------------------

n <- 300
p <- c(0.15, 0.5, 0.35)
mu1 <- 1
mu2 <- 6
mu3 <- 7
sd1 <- 3
sd2 <- 1
sd3 <- 0.5

set.seed(123)
y <- as.data.frame(r_mixNorm_2(n, mu1, mu2, mu3, sd1, sd2, sd3, p))

ggplot(data = y, mapping = aes(x = V1)) +
      geom_histogram(aes(y = after_stat(density)), 
                     fill = "gray50", color = "white", bins = 30, alpha = 0.6) +
      geom_density(color = "#1f77b4", fill = "#1f77b4", alpha = 0.3, linewidth = 1) +
      labs(title = "Simulazione 10", x = "y", y = "Densità") +
      theme_minimal(base_size = 13)

