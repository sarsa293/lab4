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
    for (i in 1:m){
      Q[k:m,i] <- Q[k:m,i]-(v*c(2*(t(v)%*%Q[k:m, i])))
    }
  }
  Q <- t(Q)
  R <- X*upper.tri(X, diag=T)
}