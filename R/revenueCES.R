#' Two-step estimation of revenue CES
#'
#' \code{revenueCES} takes production function inputs and revenue as arguments
#' and estimates a CES production function with labor-augmenting productivity
#' up to a returns-to-scale parameter, assuming only cost minimization, as in
#'  \insertCite{vandijcke2022;textual}{revenueCES}.
#'
#' Try and write
#'
#'
#'
#'
#' @references
#' \insertRef{vandijcke2022}{revenueCES} \cr
#'
#' @param R A vector of firms' revenue over time.
#' @param L A vector of firms' labor inputs over time.
#' @param M A vector of firms' material inputs over time.
#' @param P_L A vector of firm- or industry-level labor input prices.
#' @param P_M A vector of firm- or industry-level material input prices.
#' @param idvar A vector indexing which firm each row of the vectors corresponds to.
#' @param timevar A vector indexing which time period each row of the vectors corresponds to
#' @param A degree of the polynomial for the Markov process
#' @param v_bounds A Boolean indicating whether bounds on the returns-to-scale parameter
#' v should be computed. Default = FALSE.
#' @param v An optional number in \eqn{]0,\infty[} to which the returns-to-scale parameter
#' should be restricted. Constant returns to scale correspond to \code{v} = 1. Default = NULL.
#' @param mu A Boolean whether firm-level markups (with error term) should be computed and reported.
#' Requires that \code{v} is specified. Default = FALSE.
#' @param boot Number of bootstrap repetitions for calculating standard errors. Default = 100.
#' @param opt A string with the optimization algorithm to be used during the GMM estimation.
#' Options are \link[stats]{optim}, \link[Rsolnp]{solnp} and \link[DEoptim]{DEoptim}. Default = 'optim'.
#' @param theta_0 An optional named vector of starting values for the parameters,
#' where the names are 'beta_L', 'beta_M', 'sigma'. Default = NULL.
#' @param cluster An optional object of class "SOKcluster" or "cluster" for parallelization of the bootstrap.
#' See \link[paralell]{makeCluster}. Default = NULL.
#' @param control A \link[base]{list} of arguments to be passed to the optimization algorithm.
#' See the links for the optim parameter above. Default = NULL.
#' @return A member of the S3 class \code{revenueEst}, containing the following elements: \cr \cr
#' Estimates, a list with elements,
#' \item{theta}{A vector with estimated parameters.}
#' \item{theta.se}{A vector with standard errors of the estimated parameters.}
#' \item{omegal}{A vector of firm-level labor-augmenting productivity terms,
#' one for each firm in each period except the first.}
#' \item{v_bounds}{A vector with the estimated upper and lower bound on the returns to scale,
#' if specified.}
#' \item{alpha}{A vector with estimated output elasticities, if a returns-to-scale
#' restriction was specified.}
#' \item{mu}{A vector with estimated markups (with error term), if requested and
#' a returns-to-scale restriction was specified.}
#' Model, a list with elements,
#' \item{theta0}{Numeric object with the optimization starting points.}
#' \item{thetaFS}{Numeric object with the first-step parameters.}
#' \item{boot.repetitions}{The number of bootstrap repetitions used for computation of the
#' standard errors.}
#' \item{elapsed.time}{Time elapsed during the estimation.}
#' \item{opt}{String with the optimization routine used.}
#' \item{opt.outcome}{Optimization outcome.}
#' Data, a list with elements,
#' \item{R}{A vector of firms' revenue over time.}
#' \item{L}{The vector of firms' labor inputs}
#' \item{L}{The vector of firms' material inputs}
#' \item{P_L}{A vector of firm- or industry-level labor input prices.}
#' #\item{P_M}{A vector of firm- or industry-level material input prices.}
#' \item{idvar}{A vector indexing which firm each row of the vectors corresponds to.}
#' \item{timevar}{A vector indexing which time period each row of the vectors corresponds to.}
#'
#' @examples
#' Example
revenueCES <- function(R, L, M, P_L, P_M, S_L, S_M, idvar, timevar, A = 3, v_bounds = FALSE, v = NULL, mu = NULL, boot = 100, opt = 'optim',
                       theta_0 = NULL, cluster = NULL, control = NULL){

  set.seed(seed)
  Start = Sys.time() # start tracking time
  R <- checkMat(R) # change all input to matrix
  L <- checkMat(L)
  M <- checkMat(M)
  P_L <- checkMat(P_L)
  P_M <- checkMat(P_M)
  S_L <- checkMat(S_L)
  S_M <- checkMat(S_M)
  idvar <- checkMat(idvar)
  timevar <- checkMat(timevar)

  if (!all(lapply(list(R, L, M, P_L, P_M, idvar, timevar), length) == length(L))) stop("Data vectors are not of identical length") # check that all vectors have same length

  if (length(theta0) != 3){
    stop(paste0('theta0 length (', length(theta0), ') is inconsistent with the number of parameters in the revenue CES (3)'), sep = '')
  }

  lag.R = lagPanel(R, idvar = idvar, timevar = timevar)
  lag.L = lagPanel(L, idvar = idvar, timevar = timevar)
  lag.M = lagPanel(M, idvar = idvar, timevar = timevar)
  lag.P_L = lagPanel(P_L, idvar = idvar, timevar = timevar)
  lag.P_M = lagPanel(P_M, idvar = idvar, timevar = timevar)

  # set up time fixed effects
  t <- uniqueN(timevar)
  t_mat <- fastDummies::dummy_cols(timevar)
  t_mat$.data <- NULL
  t_mat <- as.matrix(t_mat)

  data <- as.matrix(model.frame(R, L, M, P_L, P_M, lag.L, lag.M, lag.P_L, lag.P_M, S_L, S_M, t_mat,
                                Z = data.frame(lag.L, lag.R, lag.M, lag.P_L, lag.P_M, t_mat)))

  # get coefficient estimates
  betas <- finalDVD(ind = TRUE, data, idvar = idvar,  opt = opt, theta0 = theta0, A = A)

  # get labor productivity estimates
  sigma <- betas$betas[1]
  beta_L <- betas$betas[2]
  beta_M <- betas$betas[3]
  omegal <- (beta_M/beta_L)^(1/sigma) * (P_L / P_M)^(1/sigma) * (L / M)^((1-sigma)/sigma)

  if (betas$opt.outcome$convergence != 0){
    warning('Second Stage convergence not achieved')
  }
  boot.indices <- block.boot.resample(idvar, boot) # generates a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    nCores = NULL
    boot.betas <- matrix(unlist(
      lapply(boot.indices, finalDVD, data = data, idvar = ivar, opt = opt,
             theta0 = theta0, A = A, boot = TRUE)), ncol = 3, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  } else {
    nCores = length(cluster)
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix( unlist( parLapply(cl = cluster, boot.indices, finalDVD, data = data, ivar = idvar,  opt = opt,
                                            theta0 = theta0, A = A, boot = TRUE)), ncol = 3, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas, 2, sd, na.rm = TRUE) # calculate standard deviations

  res.names <- c("sigma", "beta_L", "beta_M")
  names(betas$betas) <- res.names # change results' names
  names(boot.errors) <- res.names # change results' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'DVD', FSbetas = NA, boot.repetitions = boot, elapsed.time = elapsedTime, theta0 = theta0,
                          opt = opt, seed = seed, opt.outcome = betas$opt.outcome, nCores = nCores),
             Data = list(R = R, L = L, M = M, P_L = P_L, P_M = P_M, S_L = S_L, S_M = S_M, idvar = idvar, timevar = timevar),
             Estimates = list(pars = betas$betas, std.errors = boot.errors, omegal = omegal))
  return(out)
}



# function to estimate and bootstrap ACF #
finalDVD <- function(ind, data, idvar, opt, theta0, A, boot = FALSE){

  ## first estimate moment on labor-augmenting productivity

  if (sum(as.numeric(ind)) == length(ind)){ # if the ind variable is not always TRUE
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  data <- data[ind,] # change the index according to bootstrapped indices
  if (is.null(theta0)) {
    theta0 <- c(0.6) # beta_l, beta_m, sigma
  } # use reasonable guesses as starting points if user didn't specify starting


  newtime <- data[,'timevar', drop = FALSE]
  rownames(newtime) <- NULL
  J <- nrow(data[,'Z'])
  W <- matrix(1, J, J) # use identity weighting matrix

  # add time fixed effects
  beta0 <- c(theta0, rep(1, uniqueN(newtime)))

  if (opt == 'optim'){
    try.out <- try(optim(beta0, gDVD, method = "BFGS", mZ = data[, 'Z'], mW = W,
                         mL = data[, 'L'], mM = data[, 'M'], mP_L = data[, 'P_L'],
                         lag.mL = data[, 'lag.L'], lag.mM = data[, 'lag.M'],
                         lag.mP_L = data[, 'P_L'], tmat = data[, 't_mat'],
                         mA = A, control = control), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$par
      sigma <- betas[1] # drop time fixed effects, retain sigma
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA, 3, 1)
      opt.outcome <- list(convergence = 999)
      print("First-stage convergence failed")
      return(list(betas = betas, opt.outcome = opt.outcome))
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'DEoptim'){
    if (!is.null(control)) { control$trace <- FALSE } # don't print optimization output
    try.out <- try(DEoptim(gDVD, lower =  rep.int(0,length(beta0)), upper = rep.int(5,length(beta0)),
                           mZ = data[, 'Z'], mW = W,
                           mL = data[, 'L'], mM = data[, 'M'], mP_L = data[, 'P_L'],
                           lag.mL = data[, 'lag.L'], lag.mM = data[, 'lag.M'],
                           lag.mP_L = data[, 'P_L'], tmat = data[, 't_mat'],
                           mA = A, control = control), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$par
      sigma <- betas[1] # drop time fixed effects, retain sigma
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA, 3, 1)
      opt.outcome <- list(convergence = 999)
      print("First-stage convergence failed")
      return(list(betas = betas, opt.outcome = opt.outcome))
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'solnp'){
    if (!is.null(control)) { control <- list(trace = FALSE)}  # don't print optimization output
    try.out <- try(solnp(beta0, gDVD, mZ = data[, 'Z'], mW = W,
                         mL = data[, 'L'], mM = data[, 'M'], mP_L = data[, 'P_L'],
                         lag.mL = data[, 'lag.L'], lag.mM = data[, 'lag.M'],
                         lag.mP_L = data[, 'P_L'], tmat = data[, 't_mat'],
                         mA = A, control = control), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$par
      sigma <- betas[1] # drop time fixed effects, retain sigma
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA, 3, 1)
      opt.outcome <- list(convergence = 999)
      print("First-stage convergence failed")
      return(list(betas = betas, opt.outcome = opt.outcome))
    } # error handling: if the optimization fails return missing values
  }


  ## now use estimates to get beta_L and beta_M from production function
  # V = L
  pf_M <- pfDVD(mSigma = sigma, mR = data[,'R'], mL = data[, 'L'], mM = data[, 'M'], mP_L = data[, 'P_L'],
                mP_M = data[, 'P_M'], V = data[, 'M'], S_V = data[, 'S_M'])
  mdl <- lm(pf_M ~ 1)
  beta_M <- mdl$coefficients

  # V = M
  pf_L <- pfDVD(mSigma = sigma, mR = data[,'R'], mL = data[, 'L'], mM = data[, 'M'], mP_L = data[, 'P_L'],
                mP_M = data[, 'P_M'], V = data[, 'L'], S_V = data[, 'S_L'])
  mdl <- lm(pf_L ~ 1)
  beta_L <- mdl$coefficients

  betas <- c('sigma' = sigma, 'beta_L' = beta_L, 'beta_M' = beta_M)

  if (boot == FALSE){ # for the baseline estimation we want info on optimization, too
    return(list(betas = betas, opt.outcome = opt.outcome))
  } else {
    return(betas)
  }
}
# end of ACF final function #


gDVD <- function(beta, mZ, mW, mL, mM, mP_L, lag.mL, lag.mM, lag.mP_L, tmat, mA){
  # function to run the GMM estimation for ACF #

  sigma <- beta[1]
  delta <- beta[2:(2+ncol(tmat)-1)] # time fixed effects
  Omega <- delta * tmat + mP_L + (1/sigma) * (mP_L) + (1-sigma) / sigma * (mL - mM)
  Omega_lag <-  delta * tmat + mP_L + (1/sigma) * (mP_L)  + (1-sigma) (lag.mL - lag.mM) / sigma
  Omega_lag_pol <- poly(Omega_lag, degree = mA, raw=TRUE) # create polynomial in lagged omega for given degree A
  Omega_lag_pol <- cbind(1, Omega_lag_pol)
  g_b <- solve(crossprod(Omega_lag_pol)) %*% t(Omega_lag_pol) %*% Omega # OLS solution to get AR(A) coefficients for omega
  XI <- Omega - Omega_lag_pol %*% g_b # shock to productivity
  crit <- t(crossprod(mZ, XI)) %*% mW %*% (crossprod(mZ, XI))
  return(crit)

}
# end of GMM ACF #


pfDVD <- function(mSigma, mR, mL, mM, mP_L, mP_M, V, S_V) {
  # function to set up second step of estimation using revenue PF #

  omega_tilde <- (mP_L / mP_M)^(1/sigma) * (mL / mM)^((1-sigma)/sigma)
  lhs <- 1/2*(R - sigma * V + (1-sigma)/sigma * log( (omega_tilde*mL)^sigma + mM^sigma ) +
    (sigma-1)/sigma * log( (mP_L / omega_tilde)^(sigma/(sigma-1)) + mP_M^(sigma/(sigma-1))) - log(S_V))
  return(lhs)

}
