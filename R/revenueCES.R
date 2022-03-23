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
#' @param Y A vector of firms' revenue over time.
#' @param X A vector of state / dynamic inputs (e.g. capital).
#' @param V A vector of variable inputs (e.g. material inputs).
#' @param P_V A vector of firm- or industry-level input prices, ranked in the
#' same way as the corresponding inputs in \code{V}.
#' @param idvar A vector indexing which firm each row of the vectors corresponds to.
#' @param timevar A vector indexing which time period each row of the vectors corresponds to
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
#' \item{Y}{A vector of firms' revenue over time.}
#' \item{K}{The vector of firms' capital inputs}
#' \item{L}{The vector of firms' capital inputs}
#' \item{P_V}{A vector of firm- or industry-level input prices.}
#' \item{idvar}{A vector indexing which firm each row of the vectors corresponds to.}
#' \item{timevar}{A vector indexing which time period each row of the vectors corresponds to.}
#'
#' @examples
#' Example
revenueCES <- function(Y, K, L, M, P_L, P_M,idvar, timevar, v_bounds = FALSE, v = NULL, mu = NULL, boot = 100, opt = 'optim',
                       theta_0 = NULL, cluster = NULL, control = NULL){

  set.seed(seed)
  Start = Sys.time() # start tracking time
  Y <- prodest::checkM(Y) # change all input to matrix
  X <- prodest::checkM(X)
  V <- prodest::checkM(sX)
  pX <- checkM(pX)
  idvar <- checkM(idvar)
  timevar <- checkM(timevar)
  snum <- ncol(sX) # find the number of input variables
  fnum <- ncol(fX)
  if (!is.null(zX) & control == "2s") {zX <- checkM(zX); znum <- ncol(zX)} else {znum <- 0} # to determine number of variables later
  if (length(theta0) != znum + fnum + snum & !is.null(theta0)){
    stop(paste0('theta0 length (', length(theta0), ') is inconsistent with the number of parameters (', znum + fnum + snum, ')'), sep = '')
  }
  if (!is.null(zX)) {
    polyframe <- poly(fX,sX,pX,degree=G,raw=!orth) # create (orthogonal / raw) polynomial of degree G
    regvars <- cbind(fX,sX,zX,pX,polyframe) # to make sure 1st degree variables come first (lm will drop the other ones from 1st-stage reg)
  } else {
    polyframe <- poly(fX,sX,pX,degree=G,raw=!orth)
    regvars <- cbind(fX,sX,pX,polyframe) # to make sure 1st degree variables come first (lm will drop the other ones from 1st-stage reg)
  } # create (orthogonal / raw) polynomial of degree G
  if (dum) {regvars <- cbind(regvars, factor(timevar))} # add time dummies to first stage
  lag.sX = sX # generate sX lags
  for (i in 1:snum) {
    lag.sX[, i] = lagPanel(sX[, i], idvar = idvar, timevar = timevar)
  }
  lag.fX = fX # generate fX lags
  for (i in 1:fnum) {
    lag.fX[, i] = lagPanel(fX[, i], idvar = idvar, timevar = timevar)
  }
  if (control == "2s") { # generate zX lags only if we include controls in production function
    lag.zX = zX # generate fX lags
    for (i in 1:znum) {
      lag.zX[, i] = lagPanel(zX[, i], idvar = idvar, timevar = timevar)
    }
  }
  if (!is.null(zX) & control == "fs") { # generate the matrix of data for case where controls only appear in first stage, as in De Loecker and Warzynski (2012)
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX),
                                 Xt = data.frame(fX,sX), lX = data.frame(lag.fX,lag.sX),
                                 zX = data.frame(zX), regvars = regvars))
  } else if (!is.null(zX) & control == "2s") {
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX,lag.zX),
                                 Xt = data.frame(fX,sX,zX), lX = data.frame(lag.fX,lag.sX,lag.zX),
                                 zX = data.frame(zX), regvars = regvars))
  } else {
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX),
                                 Xt = data.frame(fX,sX), lX = data.frame(lag.fX,lag.sX), regvars = regvars))
  }


}
