library(FNN)

Iacf <- function(x, y, k_knn) {
  
  Iacf_affi <- function(x, y, delta = 1e-6) {
    # x <- as.matrix(x)
    x_pos <- x[y == 1,,drop = FALSE]
    x_neg <- x[y == -1,,drop = FALSE]
    n <- nrow(x)
    
    # Calculate the LS One-Class SVM distance for the positive class
    n_pos <- nrow(x_pos)
    e_pos <- rep(1,n_pos)
    c_pos <- n/2/n_pos
    
    H_pos <- rbind(c(0,e_pos),cbind(e_pos,x_pos%*%t(x_pos)+diag(1/c_pos, n_pos)))
    colnames(H_pos) <- NULL
    
    alpha_pos <- solve(H_pos, c(1, rep(0, n_pos)))
    d_pos <- abs(x_pos %*% t(x_pos) %*% alpha_pos[-1] + alpha_pos[1])
    
    # Calculate the LS One-Class SVM distance for the negative class
    n_neg <- nrow(x_neg)
    e_neg <- rep(1,n_neg)
    c_neg <- n/2/n_neg
    
    H_neg <- rbind(c(0,e_neg),cbind(e_neg,x_neg%*%t(x_neg)+diag(1/c_neg, n_neg)))
    colnames(H_neg) <- NULL
    
    alpha_neg <- solve(H_neg, c(1, rep(0, n_neg)))
    d_neg <- abs(x_neg %*% t(x_neg) %*% alpha_neg[-1] + alpha_neg[1])
    
    # Merge the distances and normalize the calculation of the affinity.
    s_affinity <- numeric(length(y))
    s_affinity[y == +1] <- 1 - (d_pos - min(d_pos)) / (max(d_pos) - min(d_pos) + delta)
    s_affinity[y == -1] <- 1 - (d_neg - min(d_neg)) / (max(d_neg) - min(d_neg) + delta)
    
    return(s_affinity)
  }
  Iacf_prob <- function(x, y, k_knn) {
    n <- nrow(x)
    n_ypos <- sum(y == 1)
    n_yneg <- sum(y == -1)
    
    # 使用 FNN 包计算 KNN
    knn_result <- get.knn(x, k = k_knn + 1)  # +1 因为包含自身点
    knn_indices <- knn_result$nn.index[, -1]  # 排除第一个（自身）
    
    # K <- x %*% t(x)
    # d <- sqrt(outer(diag(K), diag(K), "+") - 2 * K)
    # knn_indices <- matrix(0, n, k_knn)
    # for (i in 1:n) {
    #   knn_indices[i, ] <- order(d[i, ])[2:(k_knn + 1)]
    # }
    
    knn_labels <- matrix(y[knn_indices], n, k_knn)
    num_own <- rowSums(knn_labels == matrix(y, n, k_knn))
    
    m_probility <- 1- num_own / k_knn
    return(m_probility)
  }
  
  mu_vec <- Iacf_affi(x, y)
  prob_vec <- Iacf_prob(x, y, k_knn)
  
  # compute nu
  nu_vec <- (1 - mu_vec)*prob_vec
  
  # compute score
  score_vec <- (1 - nu_vec)/(2 - mu_vec - nu_vec)
  score_vec[which(nu_vec == 0)] <- mu_vec[which(nu_vec == 0)]
  score_vec[which(mu_vec <= nu_vec)] <- 0
  return(score_vec)
}
