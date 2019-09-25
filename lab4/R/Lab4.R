# 1.2.2 (*) Using the QR decomposition

linreg <- function(formula, data){
  linreg <- setRefClass("linreg", fields = list(y = "numeric", coefficients = "numeric", fitted_values = "numeric", residuals = "numeric", degrees_freedom = "numeric", variance = "numeric", qr = "matrix", t_values = "numeric"))
  ans <- linreg$new()
  X <- model.matrix(formula, data)
  qr <- X
  
  dvar <- all.vars(formula)[1]
  dvar <- data[[dvar]]
  ans$y <- dvar
  
  #QR Factorization starts here, reference https://youtu.be/d-yPM-bxREs 
  m <- dim(qr)[1]
  n <- dim(qr)[2]
  Q <- diag(m)
  for (k in 1:n){
    #Find the HH reflector
    z <- matrix(qr[k:m,k])
    v <- c(-sign(z[1]) * sqrt(sum(z^2)) - z[1], -z[2:length(z)])
    v <- v/c(sqrt(t(v) %*% v))
    v <- matrix(v, ncol=1)
    for (j in 1:n){
      qr[k:m,j] <- qr[k:m,j]-(v*c(2*(t(v)%*%qr[k:m, j])))
    }
    for (j in 1:m){
      Q[k:m,j] <- Q[k:m,j]-(v*c(2*(t(v)%*%Q[k:m, j])))
    }
  }
  
  #QR Matrix
  ans$qr <- qr
  
  #Getting Q and R
  Q <- t(Q)[,1:n]
  R <- (qr*upper.tri(qr, diag=T))[1:n,]
  #Calculating beta 
  right_hand_side <- t(Q) %*% dvar
  
  #Regression Coefficients
  beta_hat <- solve(R,right_hand_side)[,1] 
  ans$coefficients <- beta_hat
  
  #Fitted Values
  y_hat <- (X %*% beta_hat)[,1]
  ans$fitted_values <- y_hat
  
  #residuals
  res <- as.vector(dvar - X %*% beta_hat)
  ans$residuals <- res
  
  #degrees of freedom
  df <- m-n
  ans$degrees_freedom <- df
  
  #variance
  variance <- sum(res ^ 2) / df
  ans$variance <- variance
  
  #t-values
  t_values <- beta_hat / sqrt(variance)
  ans$t_values <- t_values
  
  return(ans)
  
}
