# 1.2.2 (*) Using the QR decomposition

#'@title Linear model using QR decompoition
#'
#'@description This function creates a linreg object. The object contains the qr matrix, coefficients, fitted values, residuals, degrees of freedom, variance, t_values and p_values of the linear model for \code{x} and \code{y}.
#'@param formula A formula with the shape x ~ y, where x is the dependent variable and y the independent variable.
#'@param data A dataset containing the variables of formula.
#'
#'@return A linreg Reference Class object containing the following statistics: 
#'@return \code{coefficients} a named vector of coefficients,
#'@return \code{residuals} the residuals, that is response minus predicted values
#'@return \code{predicted_values} the fitted mean values 
#'@return \code{variance} The variance
#'@return \code{t_values} t_value, 
#'@return \code{p_value} p_value
#'
#'@references More information of QR decomposition \href{https://en.wikipedia.org/wiki/QR_decomposition}{here}.
#'@references More information on Linear regression \href{https://en.wikipedia.org/wiki/Linear_regression}{here}.
#'
#'@export

# 1.2.2 (*) Using the QR decomposition

#Creating object
linreg <- setRefClass("linreg", fields = list(
  formula = "formula", 
  data = "data.frame", 
  X = "matrix",
  qr = "matrix",
  y = "numeric",
  Q = "matrix",
  R = "matrix",
  coefficients = "numeric",
  fitted_values = "numeric",
  residuals = "numeric",
  df= "numeric",
  variance = "numeric",
  t_values = "numeric"
))

#Modifying methods
linreg$methods(initialize = function(formula, data){
  X <<- model.matrix(formula, data)
  qr <<- X
  y <<- data[[all.vars(formula)[1]]]
  m <- dim(qr)[1]
  n <- dim(qr)[2]
  Q <<- diag(m)
  
  #Cycle to get matrix Q and QR
  for (k in 1:n){
    #Find the HH reflector
    z <- matrix(qr[k:m,k])
    v <- c(-sign(z[1]) * sqrt(sum(z^2)) - z[1], -z[2:length(z)])
    v <- v/c(sqrt(t(v) %*% v))
    v <- matrix(v, ncol=1)
    for (j in 1:n){
      qr[k:m,j] <<- qr[k:m,j]-(v*c(2*(t(v)%*%qr[k:m, j])))
    }
    for (j in 1:m){
      Q[k:m,j] <<- Q[k:m,j]-(v*c(2*(t(v)%*%Q[k:m, j])))
    }
  }
  
  # Q and R are:
  Q <<- t(Q)[,1:n]
  R <<- (qr*upper.tri(qr, diag=T))[1:n,]
  
  #Calculating beta 
  right_hand_side <- t(Q) %*% y
  
  #Regression Coefficients (beta hat)
  coefficients <<- solve(R,right_hand_side)[,1] 
  
  #Fitted Values (y hat)
  fitted_values <<- (X %*% coefficients)[,1]
  
  #residuals
  residuals <<- as.vector(y - X %*% coefficients)
  
  #degrees of freedom
  df <<- m-n
  
  #variance
  variance <<- sum(residuals ^ 2) / df
  
  #t-values
  t_values <<- coefficients / sqrt(variance)
})
#Modifying Print
linreg$methods(show = function(){print("Coefficients:"); print(coefficients)})
linreg$methods(print = function(){return(coefficients)})
#Methods for resid, pred, coef and summary
linreg$methods(resid = function(){return(residuals)})
linreg$methods(pred = function(){return(fitted_values)})
linreg$methods(coef = function(){return(coefficients)})
