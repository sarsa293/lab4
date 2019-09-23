# 1.2.2 (*) Using the QR decomposition
linreg <- function(formula, data){
  X <- model.matrix(formula, data)
  X <- unname(X)
  dvar <- all.vars(formula)[1]
  m <- dim(X)[1]
  n <- dim(X)[2]
  Q <- diag(m)
  for (k in 1:n){
    #Find the HH reflector
    z <- matrix(X[k:m,k])
    v <- c(-sign(z[1]) * sqrt(sum(z^2)) - z[1], -z[2:length(z)])
    v <- v/c(sqrt(t(v) %*% v))
    v <- matrix(v, ncol=1)
    for (j in 1:n){
      X[k:m,j] <- X[k:m,j]-(v*c(2*(t(v)%*%X[k:m, j])))
    }
    for (j in 1:m){
      Q[k:m,j] <- Q[k:m,j]-(v*c(2*(t(v)%*%Q[k:m, j])))
    }
  }
  #What is nargin == 2?? if nargin == 2 Q and R have to be truncated...
  Q <- t(Q)[,1:n]
  R <- (X*upper.tri(X, diag=T))[1:n,]
  #Calculating beta 
  right_hand_side <- t(Q) %*% y
  beta <- solve(R,right_hand_side)
}