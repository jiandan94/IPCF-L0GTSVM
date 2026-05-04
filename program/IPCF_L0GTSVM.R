source("program\\clipDCD.R")
source("program\\Iacf.R")

IPCF_L0GTSVM <- function(x,y,k_knn, tau,zeta,c1,c2,c3,kmax=5,verbose=FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  
  # w_til <- rnorm(2*p + 2)
  # d <- rnorm(2*p + 2)
  w_til <- rep(0, 2*p + 2)
  d <- rep(0, 2*p + 2)
  h <- rep(1, 2*p + 2)
  
  T <- 2 * (tau + 1)
  x_til <- cbind(x, 1)
  x_tilpos <- x_til[y == 1,]
  x_tilneg <- x_til[y == -1,]
  n_pos <- nrow(x_tilpos);n_neg <- nrow(x_tilneg)
  e_pos <- rep(1, n_pos);e_neg <- rep(1, n_neg);e <- c(e_pos,e_neg)
  
  B <- rbind(cbind(x_tilpos, x_tilpos), cbind(-x_tilneg, -x_tilneg))
  
  # Active_history <- list()
  for (k in 0:kmax) {
    Active <- which(sqrt(h) * abs(w_til + zeta * d) >= 
                      sort(abs(sqrt(h) * abs(w_til + zeta * d)), decreasing = TRUE)[T])
    Inactive <- setdiff(1:(2*p + 2), Active)
    # Active_history[[k+1]] <- Active
    
    w_til[Inactive] <- 0
    d[Active] <- 0
    
    Active_nowpos <- sort(Active[Active < p+1])
    Active_nowneg <- sort(Active[Active > p+1 & Active < 2*(p+1)] - (p+1))
    Active_now <- union(Active_nowpos, Active_nowneg)
    # Active_now <- intersect(Active_nowpos, Active_nowneg)
    
    iacfscore <- Iacf(x[, Active_now, drop = FALSE], y, k_knn)
    iacfscore[k == 0] <- 1
    iacfscore[iacfscore == 0] <- 1e-10
    
    S_pos <- diag(iacfscore[y == 1])
    S_neg <- diag(iacfscore[y == -1])
    S <- diag(c(diag(S_pos)/nrow(S_pos),diag(S_neg)/nrow(S_neg)))
    
    H <- (n/c1) * B[, Active] %*% t(B[, Active]) + (1/c3) * solve(S)
    a1 <- (1/n_neg) * t(e_neg) %*% S_neg %*% x_tilneg
    a2 <- (-1/n_pos) * t(e_pos) %*% S_pos %*% x_tilpos
    a <- cbind(a1,a2)
    l <- (n*c2/c1) * as.vector(B[, Active] %*% a[Active]) + rep(1, n)
    lb <- rep(0, n);ub <- rep(10^(3), n)
    alp <- clipDCD(H,l,lb,ub)
    w_til[Active] <- (n/c1) * (t(B[, Active]) %*% alp - c2 * a[Active])
    D <- diag(as.numeric(1 - y*(cbind(x_til,x_til)[, Active] %*% w_til[Active]) > 0))
    
    h <- c1/n + c3 * diag(t(B) %*% S %*% D %*% B) 
    d[Inactive] <- -(1/h[Inactive]) * 
      (c1/n * w_til[Inactive] + c2 * a[Inactive] - c3 * t(B[, Inactive]) 
       %*% S %*% D %*% (e - B[, Active] %*% w_til[Active]))
    
    if (k > 0) {
      if (length(Active) == length(Active_prev) && all(Active == Active_prev)) {
        break
      }
    }
    Active_prev <- Active
    # if (k > 0 && identical(active_history[[k]], active_history[[k + 1]]))break
    if (verbose && k %% 1 == 0) {
      cat("迭代", k, "：活跃集大小 =", length(Active), "\n")}
  }
  w_tilpos <- w_til[1:(p + 1)]; w_tilneg <- w_til[(p+2) : (2*p+2)]
  w_tilActivepos <- w_tilpos[w_tilpos != 0]
  w_tilActiveneg <- w_tilneg[w_tilneg != 0]
  Active_pos <- sort(Active[Active <= p+1])
  Active_neg <- sort(Active[Active > p+1 & Active <= 2*(p+1)] - (p+1))
  list(a = a, B = B, S = S,
       w_til = w_til, w_tilpos = w_tilpos, w_tilneg = w_tilneg,
       w_tilActive = w_til[Active],
       # w_tilInactive = w_til[Inactive],
       w_tilActivepos = w_tilActivepos, w_tilActiveneg = w_tilActiveneg,
       Active = Active, Inactive = Inactive,
       Active_pos = Active_pos,Active_neg = Active_neg,
       convergence = k < kmax,n_iter = k
  )
}