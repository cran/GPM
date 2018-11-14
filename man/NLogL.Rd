\name{NLogL}
\alias{NLogL}
\title{The Function for calculating the Negative Log-Likelehood in \code{GPM} Package}
\usage{
NLogL(Omega, X, Y, CorrType, MinEig, Fn, n, dy)
}
\arguments{
\item{Omega}{The vector storing all the scale (aka roughness) parameters of the correlation function. The length of \code{Omega} depends on the \code{CorrType}. See \code{reference 1}.}

\item{X}{Matrix containing the training (aka design or input) data points. The rows and columns of \code{X} denote individual observation settings and input dimension, respectively.}

\item{Y}{Matrix containing the output (aka response) data points. The rows and columns of \code{Y} denote individual observation responses and output dimension, respectively.}

\item{CorrType}{The correlation function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.}

\item{MinEig}{The smallest eigen value that the correlation matrix is allowed to have, which in return determines the appraopriate nugget that should be added to the correlation matrix.}

\item{Fn}{A matrix of \code{1}'s with \code{nrow(X)} rows and \code{1} column. See \code{reference 1}.}

\item{n}{Number of observations, \code{nrow(X)}.}

\item{dy}{Number of responses, \code{ncol(Y)}.}
}
\value{
nlogl The negative log-likelihood (excluding all the constant terms). See the \code{references}.
}
\description{
Calculates the negative log-likelihood (excluding all the constant terms) as described in \code{reference 1}.
}
\details{
\code{\link[GPM]{Fit}} calls this function with \emph{scaled} \code{X} and \code{Y}. That is, when the user fits a GP model by calling \code{Fit(X, Y)}, \code{X} and \code{Y} are mapped to the \code{[0, 1]} region and then passed to this function.
}
\note{
This function is \strong{NOT} exported once GPM package is loaded.
}
\examples{
# see the examples in the fitting function.
}
\references{
\enumerate{
\item Bostanabad, R., Kearney, T., Tao, S., Apley, D. W. & Chen, W. Leveraging the nugget parameter for efficient Gaussian process modeling. International Journal for Numerical Methods in Engineering, doi:10.1002/nme.5751.
\item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
}
}
\seealso{
\code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
\code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
\code{\link[GPM]{Draw}} to plot the response via the fitted model.
}
