# Quantile regression function
library(lpSolve)

my_rq <- function(y, X, tau){
  # Add intercept
  X <- cbind(1, X)
  
  # Store n and p
  n <- length(y)
  p <- ncol(X)
  
  # Objective function
  # relevant quantities: b+, b-, u, v
  f.obj <- c(rep(0, 2*p), tau * rep(1, n), (1 - tau) * rep(1, n))
  
  # Construct A matrix: [X, -X, I, -I]
  A.eq <- cbind(X, -X, diag(n), -diag(n))
  b.eq <- y
  
  # Solve the LP problem using lpSolve
  res <- lp("min", f.obj, A.eq, rep("=", n), b.eq)
  
  # beta =  b+ - b-
  beta <- res$solution[1:p] - res$solution[(p+1):(2*p)]
  names(beta) <- colnames(X)
  names(beta)[1] <- "(Intercept)"
  return(beta)
}
