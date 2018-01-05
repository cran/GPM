#' @title The Function for calculating the Negative Log-Likelehood in \code{GPM} Package
#' @description Calculates the negative log-likelihood (excluding all the constant terms) as described in \code{reference 1}.
#'
#' @param Omega The vector storing all the scale (aka roughness) parameters of the correlation function. The length of \code{Omega} depends on the \code{CorrType}. See \code{reference 1}.
#' @param X Matrix containing the training (aka design or input) data points. The rows and columns of \code{X} denote individual observation settings and input dimension, respectively.
#' @param Y Matrix containing the output (aka response) data points. The rows and columns of \code{Y} denote individual observation responses and output dimension, respectively.
#' @param CorrType The correlation function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.
#' @param MinEig The smallest eigen value that the correlation matrix is allowed to have, which in return determines the appraopriate nugget that should be added to the correlation matrix.
#' @param Fn A matrix of \code{1}'s with \code{nrow(X)} rows and \code{1} column. See \code{reference 1}.
#' @param k Number of observations, \code{nrow(X)}.
#' @param q Number of responses, \code{ncol(Y)}.
#'
#' @return nlogl The negative log-likelihood (excluding all the constant terms). See the \code{references}.
#'
#' @details \code{\link[GPM]{Fit}} calls this function with \emph{scaled} \code{X} and \code{Y}. That is, when the user fits a GP model by calling \code{Fit(X, Y)}, \code{X} and \code{Y} are mapped to the \code{[0, 1]} region and then passed to this function.
#'
#' @note This function is \strong{NOT} exported once GPM package is loaded.
#'
#' @references
#' \enumerate{
#' \item Bostanabad, R., Kearney, T., Tao, S., Apley, D. W. & Chen, W. Leveraging the nugget parameter for efficient Gaussian process modeling. International Journal for Numerical Methods in Engineering, doi:10.1002/nme.5751.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # see the examples in the fitting function.

NLogL <-  function(Omega, X, Y, CorrType, MinEig, Fn, k, q){

R = CorrMat(X, X, CorrType, Omega)
R <- (R + t(R))/2;

Raw_MinEig = sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
if (Raw_MinEig < MinEig){
  R = R + diag(x = 1, k, k)*(MinEig - Raw_MinEig)
}

L = t(chol(R))
if (CorrType == 'PE' || CorrType=='G'){
  FnTRinvFn = t(Fn)%*%solve(t(L), solve(L, Fn))
  B = t(Fn)%*%solve(t(L), solve(L, Y))/FnTRinvFn[1]
  temp = Y - Fn%*%B
  Sigma2 = t(temp)%*%solve(t(L), solve(L, temp))/k
  nlogl = k*log(det(Sigma2)) + q*2*sum(log(diag(L)))
}else{
  Alpha = t(Y)%*%solve(t(L), solve(L, Y))/k
  nlogl = 2*(log(det(Alpha)) + q*2*sum(log(diag(L)))/k)
}

return(nlogl)
}
