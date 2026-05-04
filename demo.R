rm(list = ls())

source("program\\IPCF_L0GTSVM.R")
library(MASS)

set.seed(1234)

# IPCF_L0GTSVM
n <- 100
p <- 10 # the total number of features
p_sig <- 2
p_redun <- p - p_sig
mu_pos <- c(seq(1, 2, 1), rep(0, p_redun))
mu_neg <- -c(seq(1, 2, 1), rep(0, p_redun))
cov_mat_nz <- matrix(-0.2, p_sig, p_sig)
cov_mat <- matrix(abs(rnorm(p*p, mean = 0.0001,sd =0.001)), p, p)
cov_mat <- (cov_mat + t(cov_mat))/2
cov_mat[1:p_sig, 1:p_sig] <- cov_mat_nz
diag(cov_mat) <- 1

x_pos <- mvrnorm(n/2, mu_pos, cov_mat)
x_neg <- mvrnorm(n/2, mu_neg, cov_mat)
x <- rbind(x_pos, x_neg)
y <- rep(c(1, -1), c(n/2, n/2))

fit <- IPCF_L0GTSVM(x, y, k_knn = 7, tau = 2,zeta = 0.5, c1 = 2^{-8},c2 = 256, c3 = 0.0625)

w_til <- (fit$w_tilpos + fit$w_tilneg)/2

y_pre <- sign(cbind(x,1) %*% w_til)
acc_results <- mean(y_pre == y)

cat(sprintf("Prediction accuracy: %.3f\n", acc_results))
