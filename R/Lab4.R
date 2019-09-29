#'@title A package for creating linreg RC objects.
#'@description This package creates a linreg object. The object contains the qr matrix, coefficients, fitted values, residuals, degrees of freedom, variance, t_values and p_values of the linear model for \code{x} and \code{y}. The object needs aid from ggplot2 to make the graphics inside linreg$plot(). 
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
#' @references More information of QR decomposition \href{https://en.wikipedia.org/wiki/QR_decomposition}{here}.
#'@references More information on Linear regression \href{https://en.wikipedia.org/wiki/Linear_regression}{here}.
#'@name lab4
#'@docType package
"_PACKAGE"
