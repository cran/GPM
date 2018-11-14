#' @title The Fitting Function of \code{GPM} Package
#' @description Fits a Gaussian process (GP) to a set of simulation data as described in \code{reference 1}. Both the inputs and outputs can be multi-dimensional. The outputs can be noisy (the noise variance is assumed to be homogeneous).
#' @param X Matrix containing the training (aka design or input) data points. The rows and columns of \code{X} denote individual observation settings and input dimension, respectively.
#' @param Y Matrix containing the output (aka response) data points. The rows and columns of \code{Y} denote individual observation responses and output dimension, respectively.
#' @param CorrType The type of the correlation function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details. For smooth (or analytic) functions, choose either \code{'G'} or \code{'LBG'}.
#' @param Eps A vector containing the smallest eigen value(s) that the correlation matrix is allowed to have. The elements of Eps must be in [0, 1] and sorted in a descending order.
#' @param Nopt The number of times the log-likelihood function is optimized when Eps[1] is used to constrain the smallest eigen value that the correlation matrix is allowed to have. Higher \code{Nopt} will increase fitting time as well as the chances of finding the global optimum.
#' @param TraceIt Non-negative integer. If positive, tracing information on the progress of the optimization is \strong{printed}. There are six levels of tracing (see \code{\link{optim}}) and higher values will produce more tracing information.
#' @param MaxIter Maximum number of iterations allowed for each optimization (see \code{\link{optim}}).
#' @param Seed An integer for the random number generator. Use this to make the results reproducible.
#' @param LowerBound,UpperBound To estimate the scale (aka roughness) parameters of the correlation function, the feasible range should be defined. \code{LowerBound} and \code{UpperBound} are vectors determining, resepectively, the lower and upper bounds and their length depends on the parametric form of the correlation function (see \code{reference 1} for the details).
#' @param StopFlag \code{Flag} indicating whether the optimization must be stopped if the constraint on the correlation matrix is inactive for two consecuitvie elements of \code{Eps}.
#' @param Progress \code{Flag} indicating if the fitting process should be summarized. Set it to \code{!=0} to turn it on.
#' @import lhs
#' @import randtoolbox
#' @return Model A list containing the following components:
#' \itemize{
#' \item{\code{CovFunc}} {A list containing the type and estimated parameters of the correlation function.}
#' \item{\code{Data}} {A list storing the original (but scaled) data.}
#' \item{\code{Details}} {A list of some parameters (used in prediction) as well as some values reporting the total run-time (\code{cost}), leave-one-out cross-validation error (\code{LOOCV}), the added nugget (\code{Nug_opt}) for satisfying the constraint on the smallest eigen value of the correlation matrix.}
#' \item{\code{Opt_History}} {The optimization history.}
#' \item{\code{Setting}} {The default/provided settings for running the code.}
#' }
#'
#' @references
#' \enumerate{
#' \item Bostanabad, R., Kearney, T., Tao, S., Apley, D. W. & Chen, W. Leveraging the nugget parameter for efficient Gaussian process modeling. International Journal for Numerical Methods in Engineering, doi:10.1002/nme.5751.
#' \item Plumlee, M. & Apley, D. W. (2017) Lifted Brownian kriging models. Technometrics, 59, 165-177.
#' }
#' @export
#' @seealso
#' \code{\link[stats]{optim}} for the details on \code{L-BFGS-B} algorithm used in optimization.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # 1D example: Fit a model (with default settings) and evaluate the performance
#' # by computing the root mean squared error (RMSE) in prediction.
#' library(lhs)
#' X <- 5*maximinLHS(8, 1)
#' Y <- 2*sin(X) + log(X+1)
#' M <- Fit(X, Y)
#' XF <- matrix(seq(0, 5, length.out = 50), 50, 1)
#' YF <- Predict(XF, M)
#' RMSE <- sqrt(mean((YF$YF - (2*sin(XF) + log(XF+1)))^2))
#'
#' \dontrun{
#'
#' # 1D example: Fit a model, evaluate the performance, and plot the response
#' # along with 95% prediction interval
#' X <- 10*maximinLHS(18, 1) - 5
#' Y <- X*cos(X)
#' M <- Fit(X, Y)
#' XF <- matrix(seq(-5, 5, length.out = 50), 50, 1)
#' YF <- Predict(XF, M)
#' RMSE <- sqrt(mean((YF$YF - (XF*cos(XF)))^2))
#' Draw(M, 1, res = 15)
#' 
#' # 2D example: Fit a model, evaluate the performance, and plot the response
#' # surface along with 95% prediction interval
#' X <- 2*maximinLHS(10, 2) - 1
#' Y <- X[, 1]^2 + X[, 2]^2
#' M <- Fit(X, Y, CorrType = "PE")
#' XF <- 2*maximinLHS(100, 2) - 1
#' YF <- Predict(XF, M)
#' RMSE <- sqrt(mean((YF$YF - (XF[, 1]^2 + XF[, 2]^2))^2))
#' library(lattice)
#' Draw(M, c(1, 1), res = 15, PI95=1)
#'
#' # 2D example: Plot the previous model wrt X1 in the [-2, 2]
#' # interval with X2=1
#' Draw(M, c(1, 0), LB = -2, UB = 2, res = 15, PI95=1)
#'
#' # 3D example: Compare the performance of Gaussian ("G") and lifted Browninan
#' # with Gamma=1 ("LBG")
#' X <- 2*maximinLHS(30, 3) - 1
#' Y <- cos(X[, 1]^2) + 2*sin(X[, 2]^2) + X[, 3]^2
#' M_G <- Fit(X, Y)
#' M_LBG <- Fit(X, Y, CorrType = "LBG")
#' XF <- 2*maximinLHS(50, 3) - 1
#' YF_G <- Predict(XF, M_G)
#' YF_LBG <- Predict(XF, M_LBG)
#' RMSE_G <- sqrt(mean((YF_G$YF - (cos(XF[, 1]^2) + 2*sin(XF[, 2]^2) + XF[, 3]^2))^2))
#' RMSE_LBG <- sqrt(mean((YF_LBG$YF - (cos(XF[, 1]^2) + 2*sin(XF[, 2]^2) + XF[, 3]^2))^2))
#'
#' # 3D example: Draw the response in 2D using the M_G model when X3=0
#' Draw(M_G, c(1, 1, 0), PI95 = 0, Values = 0, X1Label = 'Input 1', X2Label = 'Input 2')
#' }

Fit <-  function(X, Y, CorrType = 'G', Eps = 10^(seq(-1, -12)),Nopt = 5,TraceIt = 0,
                 MaxIter = 100, Seed = 1, LowerBound = NULL, UpperBound = NULL,
                 StopFlag = 0, Progress = 0) {

  set.seed(Seed)

  ## Check the data
  if (Progress != 0){
    cat('*** Checking the inputs ...\n')
  }
  if (missing(X) || missing(Y)){
    stop('     X and Y must be provided.')
  }
  if (!all(is.finite(X)) || !is.numeric(X)){
    stop('     All the elements of X must be finite numbers.')
  }
  if (!all(is.finite(Y)) || !is.numeric(Y)){
    stop('     All the elements of Y must be finite numbers.')
  }
  if (is.matrix(X) == FALSE) {
    X <- as.matrix(X)
  }
  if (is.matrix(Y) == FALSE) {
    Y <- as.matrix(Y)
  }
  if (nrow(X)!= nrow(Y)){
    stop('     The number of rows (i.e., observations) in X and Y should match!')
  }
  if (is.character(CorrType)){
    CorrType <- toupper(CorrType)
  } else {
    stop('     CorrType should be a character (options are: "G", "PE", "LB", and "LBG")')
  }
  if (CorrType!= 'G' && CorrType!= 'PE' && CorrType!= 'LBG' && CorrType!='LB'){
    stop('     The type of the correlation function is not determined correctly. Supported functions types are "G", "PE", "LB", and "LBG".')
  }
  if (any(Eps < (10^(round(log10(.Machine$double.eps))+3)))){
    stop(paste('     Increase the smallest member of Eps. The minimum allowable is ', toString(10^(round(log10(.Machine$double.eps))+3))))
  }
  if (any(diff(Eps) > 0)){
    stop('The elements of Eps should be in a descending order.')
  }

  n <- nrow(X); dx <- ncol(X); dy <- ncol(Y);


  Omega_min_default <- -2
  Omega_max_default <- 2
  Zeta_min_default <- 1e-6
  Zeta_max_default <- 1 - 1e-4
  Gamma_min_default <- 1e-6
  Gamma_max_default <- 1

  if (is.null(LowerBound)){
    LB_Def <- 1
    if (CorrType == 'G') {LowerBound = rep(x = Omega_min_default, times = dx)}
    else if (CorrType == 'PE') {LowerBound = c(rep(x = Omega_min_default, times = dx), 1)}
    else if (CorrType == 'LBG') {LowerBound = c(rep(x = Omega_min_default, times = dx), Zeta_min_default)}
    else if (CorrType == 'LB') {LowerBound = c(rep(x = Omega_min_default, times = dx), Zeta_min_default, Gamma_min_default)}
  } else {
    LB_Def <- 0
    if (!is.vector(LowerBound)) LowerBound = as.vector(LowerBound)
    if (CorrType == 'G' && length(LowerBound) != dx){
      stop(paste('     The size of the provided LowerBound is incorrect! It should be', toString(dx)))}
    else if (CorrType == 'PE' && length(LowerBound) != (1 + dx)){
      stop(paste('     The size of the provided LowerBound is incorrect! It should be', toString(dx + 1)))}
    else if (CorrType == 'LBG' && length(LowerBound) != (1 + dx)){
      stop(paste('     The size of the provided LowerBound is incorrect! It should be', toString(dx + 1)))}
    else if (CorrType == 'LB' && length(LowerBound) != (2 + dx)){
      stop(paste('     The size of the provided LowerBound is incorrect! It should be', toString(dx + 2)))}
  }

  if (is.null(UpperBound)){
    UB_Def <- 1
    if (CorrType == 'G') {UpperBound = rep(x = Omega_max_default, times = dx)}
    else if (CorrType == 'PE') {UpperBound = c(rep(x = Omega_max_default, times = dx), 2)}
    else if (CorrType == 'LBG') {UpperBound = c(rep(x = Omega_max_default, times = dx), Zeta_max_default)}
    else if (CorrType == 'LB') {UpperBound = c(rep(x = Omega_max_default, times = dx), Zeta_max_default, Gamma_max_default)}
  } else {
    UB_Def <- 0
    if (!is.vector(UpperBound)) UpperBound = as.vector(UpperBound)
    if (CorrType == 'G' && length(UpperBound) != dx){
      stop(paste('     The size of the provided UpperBound is incorrect! It should be', toString(dx)))}
    else if (CorrType == 'PE' && length(UpperBound) != (1 + dx)){
      stop(paste('     The size of the provided UpperBound is incorrect! It should be', toString(dx + 1)))}
    else if (CorrType == 'LBG' && length(UpperBound) != (1 + dx)){
      stop(paste('     The size of the provided UpperBound is incorrect! It should be', toString(dx + 1)))}
    else if (CorrType == 'LB' && length(UpperBound) != (2 + dx)){
      stop(paste('     The size of the provided UpperBound is incorrect! It should be', toString(dx + 2)))}
  }
  if (CorrType == 'PE' && (LowerBound[dx + 1] < 1 || UpperBound[dx+ 1] > 2)){
    stop(paste('     The exponent in PE should be between ', toString(c(1, 2)), ' (inclusive).'))
  } else if (CorrType == 'LBG' && (LowerBound[dx + 1] <= 0 || UpperBound[dx+ 1] > Zeta_max_default)){
    stop(paste('     Beta in LBG should be between ', toString(c(0, 1)), ' (exclusive).'))
  } else if (CorrType == 'LB' && (LowerBound[dx+1]<=0 || UpperBound[dx+1]>Zeta_max_default || 
                  LowerBound[dx+2]<Gamma_min_default || UpperBound[dx+2]>Gamma_max_default)){
    stop('     The provided range for Beta and/or Gamma is not acceptable.')
  }


  Setting <- list(CorrType = CorrType, Eps = Eps, MaxIter = MaxIter,
                 Seed = Seed, StopFlag = StopFlag, Nopt = Nopt,
                 LowerBound = LowerBound, UpperBound = UpperBound, TraceIt = TraceIt)

  ## Normalize the data
  N_hyperC <- length(LowerBound)
  Xmin <- apply(X, 2, min)
  Xmax <- apply(X, 2, max)

  XN <- t((t(X)-Xmin)/(Xmax-Xmin))
  Ymin <- apply(Y, 2, min)
  Yrange <- apply(Y, 2, max) - Ymin
  YN <- t((t(Y)-Ymin)/Yrange)
  if (CorrType == 'LBG' || CorrType == 'LB'){
    XN0 <- XN[1, ]
    YN0 <- YN[1, ]
    XN <- t(t(XN) - XN0)[2:n, ]
    YN <- t(t(YN) - YN0)[2:n, , drop = FALSE]
    n <- n - 1
  }
  Fn <- matrix(data = 1, nrow = n, ncol = 1)

  ptm <- proc.time()

  A <- maximinLHS(Nopt, N_hyperC)
  #A <- sobol(Nopt, dim = N_hyperC, scrambling = 3, seed = Seed)
  Tries <- length(Eps)
  N_unique_local_opt <- matrix(0, Tries, 1)
  OptimHist <- NULL
  OptimHist$Var_int[[1]] <- t(t(A)*(UpperBound - LowerBound) + LowerBound)
  CTRL <- c(trace = TraceIt, maxit = MaxIter,  REPORT = 1, lmm = 15, factr = 1e6, pgtol = 1e-8)

  for (i in 1:Tries){
    Var_int <- OptimHist$Var_int[[i]]
    Eps_i <- Eps[i]
    Nopt_i <- nrow(Var_int)
    Var <- matrix(0, Nopt_i, N_hyperC)
    Obj_Fun <- matrix(0, Nopt_i, 1)
    Exit_Flag <- matrix(0, Nopt_i, 1)
    if (LB_Def == 1 && i > 1){
      LowerBound[1: dx] <- -10
    }
    if (UB_Def == 1 && i > 1){
      UpperBound[1: dx] <- 4
    }
    for (j in 1: Nopt_i){
      temp <- stats::optim(Var_int[j, ], NLogL, gr = NULL, XN, YN, CorrType, Eps_i, Fn, n, dy,
                    method = 'L-BFGS-B', lower = LowerBound, upper = UpperBound, control = CTRL)
      Var[j, ] <- temp$par
      Obj_Fun[j] <- temp$value
      Exit_Flag[j] <- temp$convergence
    }
    
    ID = sort(Obj_Fun, index.return=TRUE)$ix
    Obj_Fun <- as.matrix(Obj_Fun[ID, drop = FALSE])
    Exit_Flag <- as.matrix(Exit_Flag[ID, drop = FALSE])
    Var <- as.matrix(Var[ID, , drop = FALSE])

    if (Progress != 0){
      Table <- cbind(apply(round(Var, 3), 1, toString), round(Obj_Fun, 2))
      colnames(Table)<- c('Estimated Hyperparameters', ' NLogL')
      Table <- as.table(Table)
      rownames(Table) <- 1:Nopt_i
      cat(sprintf('\n\t\tOptimization for Eps = %.1e\n', Eps[i]))
      print(Table, row.names = FALSE)
    }

    rownames(Obj_Fun) <- 1:Nopt_i
    Obj_Fun_unique  <-  as.matrix(unique(round(Obj_Fun, 2)))
    Exit_unique <- as.matrix(as.data.frame(Exit_Flag)[c(rownames(Obj_Fun_unique)), ])
    Var_unique <- as.matrix(as.data.frame(Var)[c(rownames(Obj_Fun_unique)), ])

    if (Progress != 0 && Nopt_i > 1){
      Table <- cbind(apply(round(Var_unique, 3), 1, toString), round(Obj_Fun_unique, 2))
      colnames(Table)<- c('Estimated Hyperparameters', ' NLogL')
      Table <- as.table(Table)
      rownames(Table) <- 1:nrow(Obj_Fun_unique)
      cat('\n\t\tThe unique local optima are:\n')
      print(Table, row.names = FALSE)
    }

    R <- CorrMat_Sym(XN, CorrType, Var_unique[1, ])
    Raw_MinEig <- Eigen(R)[1] #sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
    if (Raw_MinEig < Eps[i]){
      Nug_opt <- Eps[i] - Raw_MinEig;
      R <- R + diag(x = 1, n, n)*(Eps[i] - Raw_MinEig)
    } else {
      Nug_opt <- 0
    }

    OptimHist$Nug_opt[[i]] <- Nug_opt
    OptimHist$Raw_MinEig[[i]] <-  Raw_MinEig
    OptimHist$Var[[i]] <-  Var
    OptimHist$Var_unique[[i]] <-  Var_unique
    OptimHist$Obj_Fun[[i]] <-  Obj_Fun
    OptimHist$Obj_Fun_unique[[i]] <- Obj_Fun_unique
    OptimHist$Obj_Fun_best[[i]] <-  Obj_Fun[1]
    OptimHist$Exit_unique[[i]] <- Exit_unique

    if (i > 1 && (StopFlag != 0) && (OptimHist$Nug_opt[[i]] == OptimHist$Nug_opt[[i-1]])){
        break
    }
    if (i < Tries){
      OptimHist$Var_int[[i + 1]]  <-  Var_unique
    }
  }

  Cost <- proc.time()[3] - ptm[3]
  # Organize the optimization results
  #cat('*** Summary of the optimization results (sorted and rounded):\n')
  ID <- which.min(OptimHist$Obj_Fun_best)
  w <- as.matrix(OptimHist$Var[[ID]][1, ])

  ## Post-processing
  R <- CorrMat_Sym(XN, CorrType, w)
  Raw_MinEig <- Eigen(R) #sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
  if (Raw_MinEig < Eps[ID]){
    Nug_opt <- Eps[ID] - Raw_MinEig;
    if (Progress != 0){
      cat('\n*** The best model is fitted via Eps = ', toString(Eps[ID]), 
          '. A nugget of ', Nug_opt, ' has been used.')
    }
    R <- R + diag(x = 1, n, n)*(Eps[ID] - Raw_MinEig)
  } else {
    Nug_opt <- 0
  }

  L <- LowerChol(R)
  Rinv_YN <- CppSolve(t(L), CppSolve(L, YN))
  if (CorrType == 'PE' || CorrType=='G'){
    Rinv_Fn <- CppSolve(t(L), CppSolve(L, Fn))
    FnTRinvFn <- t(Fn)%*%Rinv_Fn
    B <- t(Fn)%*%Rinv_YN/FnTRinvFn[1]
    temp <- YN - Fn%*%B
    Sigma2 <- t(temp)%*%CppSolve(t(L), CppSolve(L, temp))/n
    if (CorrType=='G'){
      Theta <- w; Power <- 2;
    } else {
      Theta <- w[1:dx]; Power <- w[dx+1]
      if (Power > 1.999){
        Power <- 2
      }
    }
    Parameters <- list('FnTRinvFn' = FnTRinvFn[1], 'B' = B, 'Sigma2' = Sigma2, 'Rinv_Fn' = Rinv_Fn,
                      'Rinv_YN' = Rinv_YN, 'Theta' = Theta, 'Power' = Power)
  } else {
    Alpha <- t(YN)%*%Rinv_YN/n
    if (CorrType == 'LBG'){
      A <- w[1:dx]; Beta <- w[dx+1]; Gamma <- 1
    } else {
      A <- w[1:dx]; Beta <- w[dx+1]; Gamma <- w[dx+2]
      if (Gamma>0.999) Gamma <- 1
    }
    Parameters <- list('XN0' = XN0, 'YN0' = YN0, 'Alpha' = Alpha, 'Rinv_YN' = Rinv_YN,
                      'A' = A, 'Beta' = Beta, 'Gamma' = Gamma)
  }
  ## Save the results
  Model <- NULL
  Model$CovFunc <- list('CorrType' = CorrType, 'Parameters' = Parameters, 
                        'LowerBound' = LowerBound, 'UpperBound' = UpperBound)
  Model$Data <- list('XN' = XN, 'n' = n, 'Xmin' = Xmin, 'Xmax' = Xmax,  'YN' = YN, 
                     'dy' = dy, 'Yrange' = Yrange, 'Ymin' = Ymin)
  Model$Details <- list('Fn' = Fn, 'L' = L, 'Raw_MinEig' = Raw_MinEig, 'Nug_opt' = Nug_opt, 
                        'Cost' = Cost, 'MinNLogL' = OptimHist$Obj_Fun_best[ID])
  Model$OptimHist <- OptimHist
  Model$Setting <- Setting
  class(Model) <- 'GPM'

  return(Model)
}
