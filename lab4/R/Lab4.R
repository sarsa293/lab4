# 1.2.2 (*) Using the QR decomposition
linreg <- function(formula, data){
  X <- model.matrix(formula, data)
  X <- unname(X)
  x <- X
  dvar <- all.vars(formula)[1]
  dvar <- data[[dvar]]
  dvar <- matrix(dvar, ncol=1)
  
  #QR Factorization starts here, reference https://youtu.be/d-yPM-bxREs 
  m <- dim(x)[1]
  n <- dim(x)[2]
  Q <- diag(m)
  for (k in 1:n){
    #Find the HH reflector
    z <- matrix(x[k:m,k])
    v <- c(-sign(z[1]) * sqrt(sum(z^2)) - z[1], -z[2:length(z)])
    v <- v/c(sqrt(t(v) %*% v))
    v <- matrix(v, ncol=1)
    for (j in 1:n){
      x[k:m,j] <- x[k:m,j]-(v*c(2*(t(v)%*%x[k:m, j])))
    }
    for (j in 1:m){
      Q[k:m,j] <- Q[k:m,j]-(v*c(2*(t(v)%*%Q[k:m, j])))
    }
  }
  #What is nargin == 2?? if nargin == 2 Q and R have to be truncated...
  Q <- t(Q)[,1:n]
  R <- (x*upper.tri(x, diag=T))[1:n,]
  #Calculating beta 
  right_hand_side <- t(Q) %*% dvar
  beta_hat <- solve(R,right_hand_side) #Coefficients
  #residuals
  res <- dvar - X %*% beta
  #degrees of freedom
  df <- m-n
  #variance
  variance <- sum(res ^ 2) / df
  return(beta)
  #beta_hat <- solve(t(X) %*% x) %*% t(X) %*% y
  
  #linreg <- setRefClass("linreg", fields = list( formula = "formula", data = "data.frame", beta = "numeric"))
  
}
