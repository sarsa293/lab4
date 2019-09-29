# 1.2.2 (*) Using the QR decomposition

#'
#'@title lab4: A package for creating linreg RC objects.
#'
#'@description This package creates a linreg object. The object contains the qr matrix, coefficients, fitted values, residuals, degrees of freedom, variance, t_values and p_values of the linear model for \code{x} and \code{y}.
#'
#' @field formula A formula object like x ~ y.
#' @field data A data frame which contains the formula variables.
#' @field X The matrix created from formula.
#' @field qr The QR matrix
#' @field y A vector with the independent variable values
#' @field Q The matrix Q
#' @field R The matrix R
#' @field coefficients A vector with the beta coefficients
#' @field fitted_values A vector with the predicted values of y.
#' @field residuals A vector with the residuals (y - fitted_values).
#' @field df The degrees of freedom.
#' @field Bvar A vector containing the variance of the coefficients.
#' @field resstdev A vector with the standard deviation of the residuals
#' @field t_values A vector with the t_values of the coefficients.
#' @field p_values A vector with the p_values of the coefficients.
#' 
#'@references More information of QR decomposition \href{https://en.wikipedia.org/wiki/QR_decomposition}{here}.
#'@references More information on Linear regression \href{https://en.wikipedia.org/wiki/Linear_regression}{here}.
#'@name lab4
#'@docType package
#'@export
NULL

# 1.2.2 (*) Using the QR decomposition

library(ggplot2)

#' A Reference Class to represent a linear model.
#'
#' @field formula A formula object like x ~ y.
#' @field data A data frame which contains the formula variables.
#' @field X The matrix created from formula.
#' @field qr The QR matrix
#' @field y A vector with the independent variable values
#' @field Q The matrix Q
#' @field R The matrix R
#' @field coefficients A vector with the beta coefficients
#' @field fitted_values A vector with the predicted values of y.
#' @field residuals A vector with the residuals (y - fitted_values).
#' @field df The degrees of freedom.
#' @field Bvar A vector containing the variance of the coefficients.
#' @field resstdev A vector with the standard deviation of the residuals
#' @field t_values A vector with the t_values of the coefficients.
#' @field p_values A vector with the p_values of the coefficients.
#' @export

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
  Bvar = "numeric",
  resstdev = "numeric",
  t_values = "numeric",
  p_values = "numeric",
  data_print = "character",
  starv = "character"
))

#Modifying methods
linreg$methods(initialize = function(formula, data){
  X <<- model.matrix(formula, data)
  qr <<- X
  y <<- data[[all.vars(formula)[1]]]
  data_print <<- tail(all.vars(formula))
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
  
  #residualvariance and standard deviation
  resvariance <- as.numeric((t(residuals) %*% residuals) / df)
  resstdev <<- sqrt(resvariance)
  
  #variance of beta coefficients
  Bvar <<- diag(resvariance * solve(t(X) %*% X))
  
  #t-values
  t_values <<- coefficients / sqrt(Bvar)
  
  #pvalues
  p_values <<- 2 * pt(abs(t_values), df, lower.tail = FALSE)
  
  for(i in 1:length(coefficients)){
    starv[i] <<- "***"
  }
  
  #strings formula and data
  formula <<- formula
  data_print <<- deparse(substitute(data))
})
#Modifying Print
linreg$methods(show = function(){print("Coefficients:"); print(coefficients)})
linreg$methods(print = function() {
  cat(paste("linreg(formula = ", format(formula),", data = ", format(data_print),")\n\nCoefficients:\n", sep = ""))
  print.table(coefficients)
})
#Methods for resid, pred, coef
linreg$methods(resid = function(){return(residuals)})
linreg$methods(pred = function(){return(fitted_values)})
linreg$methods(coef = function(){return(coefficients)})
#Summary
linreg$methods(summary = function (){
  summarytable <- data.frame("Coefficients" = coefficients, "Standard error" = sqrt(Bvar), "T Values" = t_values, "P Values"= p_values, "stars" = starv)
  names(summarytable)[5] <- " "
  print.data.frame(summarytable)
  cat(paste("\nResidual standard error: ", resstdev, " on ", df, " degrees of freedom", sep = ""))
})
linreg$methods(plot = function(){
  a <- ggplot(data, aes(x = fitted_values , y = residuals))+ 
    geom_point(shape = 1) + 
    ggtitle("Residuals vs Fitted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    labs( x = "Fitted values \n linreg(Petal.Length ~ Species)" , y = "Residuals") +
    stat_summary(fun.y = median, geom = "smooth", se=TRUE, size = 1.2)
   
  b <-  ggplot(data, aes(x = fitted_values , y = sqrt(abs(residuals / sqrt(Bvar))))) + 
    geom_point(shape = 1) +
    ggtitle("Scale-Location") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs( x = "Fitted values \n linreg(Petal.Length ~ Species)" , y = "Standardized Residuals") +
    stat_summary(fun.y = median, geom = "smooth", se=TRUE, size = 1.2)
  
  list(a , b)          
  
  })

#linreg_mod <- linreg$new(Petal.Length ~ Species, data=iris)
#linreg_mod$plot()

#linreg_mod <- linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)
#linreg_mod$summary()
