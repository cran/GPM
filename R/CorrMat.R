#' @title The Function for Constructing the Correlation Matrix in \code{GPM} Package
#' @description Builds the correlation matrix given two datasets, and the type and parameters of the correlation function.
#' @param X1,X2 Matrices containing the data points. The rows and columns of both \code{X1} and \code{X2} denote individual observation settings and dimension, respectively.
#' @param CorrType The correlation function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.
#' @param Omega The vector storing all the scale (aka roughness) parameters of the correlation function. The length of \code{Omega} depends on the \code{CorrType}. See \code{reference 1}.
#'
#' @return R The Correlation matrix with size \code{nrow(X1)}-by-\code{nrow(X2)}. See \href{https://en.wikipedia.org/wiki/Correlation_matrix}{here}.
#'
#' @note This function is \strong{NOT} exported once the GPM package is loaded.
#' @references
#' \enumerate{
#' \item Bostanabad, R., Kearney, T., Tao, S. Y., Apley, D. W. & Chen, W. (2018) Leveraging the nugget parameter for efficient Gaussian process modeling. International Journal for Numerical Methods in Engineering, 114, 501-516.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # see the examples in \code{\link[GPM]{Fit}}

CorrMat <-  function(X1, X2, CorrType, Omega){
k = nrow(X1)
p = ncol(X1)
m = nrow(X2)
R = matrix(0, k, m)
Omega = as.vector(Omega)
if (CorrType == 'G'){
  if (p != length(Omega)){
    stop(paste('In a Gaussian Correlation function, there should be', toString(p), 'parameters!'))
  }
  if (k >= m){
    for (i in 1: m) {
      R[, i] = colSums((t(X1) - X2[i, ])^2*(10^Omega))
    }
    R = exp(-R)
  } else{
    for (i in 1: k) {
      R[i , ] = colSums((t(X2) - X1[i, ])^2*(10^Omega))
    }
    R = exp(-R)
  }
} else if (CorrType == 'PE'){
  if (p != length(Omega) - 1){
    stop(paste('In a PE Correlation function, there should be', toString(p+1), 'parameters!'))
  }
  if (k >= m){
    for (i in 1: m) {
      R[ , i] = exp(-colSums(abs((t(X1) - X2[i, ]))^Omega[p+1]*(10^Omega[-(p+1)])))
    }
  } else{
    for (i in 1: k) {
      R[i , ] = exp(-colSums(abs((t(X2) - X1[i, ]))^Omega[p+1]*(10^Omega[-(p+1)])))
    }
  }
} else if (CorrType == 'LBG'){
  if (p != length(Omega) - 1){
    stop(paste('In an LBG Correlation function, there should be', toString(p+1), 'parameters!'))
  }
  A = 10^Omega[1:p]
  Beta = Omega[p+1]
  if (k >= m){
    for (i in 1: m) {
      R[, i] = (1 + colSums(t(X1**2)*A))**Beta + (1 + sum(X2[i, ]**2*A))**Beta -
        (1 + colSums((t(X1) - X2[i, ])**2*A))**Beta - 1
    }
  } else{
    for (i in 1: k) {
      R[i, ] = (1 + colSums(t(X2**2)*A))**Beta + (1 + sum(X1[i, ]**2*A))**Beta -
        (1 + colSums((t(X2) - X1[i, ])**2*A))**Beta - 1
    }
  }
} else if (CorrType == 'LB'){
  if (p != length(Omega) - 2){
    stop(paste('In an LB Correlation function, there should be', toString(p+2), 'parameters!'))
  }
  A = 10^Omega[1:p]
  Beta = Omega[p+1]
  Gamma = Omega[p+2]
  if (k >= m){
    for (i in 1: m) {
      R[, i] = (1 + colSums(t(X1**2)*A)**Gamma)**Beta + (1 + sum(X2[i, ]**2*A)**Gamma)**Beta -
        (1 + colSums((t(X1) - X2[i, ])**2*A)**Gamma)**Beta - 1
    }
  } else{
    for (i in 1: k) {
      R[i, ] = (1 + colSums(t(X2**2)*A)**Gamma)**Beta + (1 + sum(X1[i, ]**2*A)**Gamma)**Beta -
        (1 + colSums((t(X2) - X1[i, ])**2*A)**Gamma)**Beta - 1
    }
  }
} else {
  stop('The type of the Correlation is not supported!')
}
R[R < (.Machine$double.eps)] <- 0

return(R)
}
