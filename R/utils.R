
simAR1 <- function(rho, sd, X, t, b0 = 0) {
  #' simulate an AR1 process forward for a vector
  #' @param rho coefficient vector
  #' @param X vector to be simulated forward t periods
  #' @param t number of periods to simulate forward

  Xmat <- as.matrix(X)
  Xlast <- X
  r0 <- X*(1-rho) # take initial value as mean
  for (step in 1:(t-1)) {
    Xnext <- r0 + rho*X + rnorm(nrow(as.matrix(Xlast)), mean = 0, sd = sd)
    Xmat <- cbind(Xmat, as.matrix(Xnext))
    Xlast <- Xnext
  }

  return(Xmat)

}


checkMat <- function(input){ # , inputname = NA
  if (is.numeric(input)) {
    stop(paste0("'", deparse(substitute(vec)), "' not numeric data. Please pass as numeric vector or matrix."))
    }
  if (!is.matrix(input)) {
    out <- as.matrix(input)
  } else{
    out <- input
  }
  colnames(out) <- NULL
  return(out)
}


