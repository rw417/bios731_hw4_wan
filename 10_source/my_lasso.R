# My LASSO function using LowRankQP
library(LowRankQP)

my_lasso <- function(y, X, lambda){
  # Add intercept
  X <- cbind(1, X)
  X <- as.matrix(X)
  p <- ncol(X)
  output_names <- colnames(X)

  # Compute matrices
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  
  # Quadratic term: (beta^+ - beta^-)^T XtX (beta^+ - beta^-)
  # Let theta = [beta^+, beta^-]^T
  # then: beta = (beta^+ - beta^-) = (I -I) * theta
  # Dimensions:
  # beta = (beta^+ - beta^-): p x 1
  # theta: 2p x 1
  # (I -I): p x 2p
  # XtX: p x p
  # (beta^+ - beta^-)^T XtX (beta^+ - beta^-)
  # = theta^T (I -I)^T XtX (I -I) theta
  II <- cbind(diag(p), -diag(p))
  Vmat <- t(II) %*% XtX %*% II
  
  # Linear term: -2 y^T X (beta^+ - beta^-)
  # = -2 y^T X (I -I) theta
  # dvec = -(y^T X (I -I))^T = (I -I)^T X^T y
  dvec <- -t(II) %*% Xty
  
  # Inequality constraint: sum(beta^+ + beta^-) <= lambda
  # So, A^T theta <= b
  # Intercept term is not constrained, so we also need to set the first term
  # and the p+1 th term of A to 0 instead of 1.
  A <- t(matrix(rep(c(0, rep(1, p - 1)), 2)))
  b <- lambda
  
  # Non-negativity constraints: theta >= 0 → enforced via uvec
  uvec <- rep(10000, 2 * p)
  
  # Solve QP
  result <- LowRankQP(Vmat, dvec, A, b, uvec, method="LU", niter=1000)
  
  # Recover beta = beta^+ - beta^-
  beta_plus  <- result$alpha[1:p]
  beta_minus <- result$alpha[(p+1):(2*p)]
  beta_hat <- beta_plus - beta_minus
  names(beta_hat) <- output_names
  
  # Set values <= 1e-6 to 0
  beta_hat[abs(beta_hat) <= 1e-6] <- 0
  return(beta_hat)
}

# Cross validation to select the best constraint for my LASSO function
# CV method = MSE
cv.my_lasso <- function(y, X, lambda_seq=NULL, k = 10) {
  # Create a sequence of lambda values
  if(is.null(lambda_seq)){
    lambda_seq <- seq(0, 150, length=100)
  }
  
  n <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  cv_errors <- numeric(length(lambda_seq))
  
  # Store the coefficients
  beta_hats <- matrix(NA, ncol=ncol(X)+1, nrow=length(lambda_seq))
  colnames(beta_hats) <- c("Intercept", colnames(X))
  
  # Loop over lambda values
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    mse_fold <- numeric(k)
    
    for (fold in 1:k) {
      test_idx <- which(folds == fold)
      train_idx <- setdiff(1:n, test_idx)
      
      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test  <- cbind(1, X[test_idx, , drop = FALSE])
      y_test  <- y[test_idx]
      
      cv_lasso_coef <- my_lasso(y_train, X_train, lambda)
      y_pred <- X_test %*% matrix(cv_lasso_coef)
      mse_fold[fold] <- mean((y_test - y_pred)^2)
    }
    
    cv_errors[i] <- mean(mse_fold)
    beta_hats[i,] <- my_lasso(y, X, lambda)
  }
  
  best_idx <- which.min(cv_errors)
  list(
    lambda_seq = lambda_seq,
    cv_errors = cv_errors,
    best_lambda = lambda_seq[best_idx],
    beta_hats = beta_hats,
    beta_best = beta_hats[best_idx,],
    best_idx = best_idx
  )
}

# 
# library(glmnet)
# lasso_glmnet = cv.glmnet(
#   x=as.matrix(cannabis[,-c(1,2)]),
#   y=log(cannabis$t_mmr1+1e-6),
#   family="gaussian",
#   alpha=1
#   )
# lasso_glmnet$lambda.min
# coef(lasso_glmnet, s=lasso_glmnet$lambda.min)
# # Corresponding coefficients (excluding intercept)
# beta <- coef(lasso_glmnet, s = "lambda.min")[-1]  # remove intercept
# # Compute L1 norm
# sum(abs(beta))
# 
# # Compute L1 norm for all penalty lambdas
# fit <- lasso_glmnet$glmnet.fit
# l1_norms <- apply(abs(as.matrix(fit$beta)), 2, sum)
# l1_norms
# 
# l1_norms <- apply(abs(as.matrix(lasso_glmnet$beta)), 2, sum)
# desired_t <- 5  # example L1-norm constraint
# lambda_index <- which.min(abs(l1_norms - desired_t))
# lasso_glmnet$beta[, lambda_index]; l1_norms[lambda_index]
# 
# cv.my_lasso(y=log(cannabis$t_mmr1+1e-6), X=as.matrix(cannabis[,-c(1,2)]),
#             lambda_seq=l1_norms, k=10)
