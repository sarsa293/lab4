# 1.2.2 (*) Using the QR decomposition
linreg <- function(formula, data){
  X <- model.matrix(formula, data)
  dvar <- all.vars(formula)[1]
  m <- dim(X)[1]
  n <- dim(X)[2]
}