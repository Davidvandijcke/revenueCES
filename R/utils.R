
simAR1 <- function(rho, sd, X, t, b0 = 0) {
  #' simulate an AR1 process forward for a vector
  #' @param rho coefficient vector
  #' @param X vector to be simulated forward t periods
  #' @param t number of periods to simulate forward

  Xmat <- as.matrix(X)
  Xlast <- X
  for (step in 1:(t-1)) {
    Xnext <- pracma::repmat(b0, length(X), 1) + rho*X + rnorm(nrow(as.matrix(Xlast)), mean = 0, sd = sd)
    Xmat <- cbind(Xmat, as.matrix(Xnext))
    Xlast <- Xnext
  }

  return(Xmat)

}
