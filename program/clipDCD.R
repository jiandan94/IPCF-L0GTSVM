
# 定义clipDCD函数
clipDCD <- function(H, l, lb, ub, eps=0.0001) {
  # H is a n*n matrix
  # l is a n*1 matrix
  # lb is a length-n vector of lower bound
  # ub is a length-n vector of upper bound
  
  # initialization
  u <- (lb + ub) / 2
  # e <- matrix(1, nrow = nrow(H), ncol = 1)
  stop_cond <- 1
  while (stop_cond >= eps) {
    # compute lambda candidate
    lam_cand <- (c(l) - t(u) %*% H) / diag(H)
    # find index set A
    A1_1 <- which(u > lb)
    A1_2 <- which(lam_cand < 0)
    A1 <- intersect(A1_1, A1_2)
    
    A2_1 <- which(u < ub)
    A2_2 <- which(lam_cand > 0)
    A2 <- intersect(A2_1, A2_2)
    A <- c(A1, A2)
    
    # 检查索引集A是否为空
    if (length(A) == 0) {
      break
    }
    
    # find L 
    decent_value <- (c(l) - t(u) %*% H)^2 / diag(H)
    decent_value_cand <- decent_value[A]
    maxind_decent_value <- which.max(decent_value_cand)
    L <- A[maxind_decent_value]
    
    # compute lam
    lam <- (l[L] - c(t(u) %*% H[, L])) / diag(H)[L]
    # update u
    u[L] <- u[L] + lam
    u[L] <- max(lb[L], min(u[L], ub[L]))
    
    # compute stop_cond
    stop_cond <- (l[L] - c(t(u) %*% H[, L]))^2 /(2*diag(H)[L])
  }
  return(u)
}