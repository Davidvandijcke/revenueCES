#as.matrix(unlist(lapply(test, function(x) matrix(x, 5,1))))

simulation <- function() {

  # assign calibrated parameters to global scope
  beta_L <<- 0.4
  beta_M <<- 0.4
  v <<- 1 # constant returns to scale
  sigma_CES <<- 0.5 # Raval (2019)
  Sigma <<- 1 # as in De Ridder et. al.
  epsi <<- 10 # as in De Ridder et. al.
  N <<- 180 # as in De Ridder et. al.
  N_h <<- 8 # as in De Ridder et. al.
  N_f <<- N*N_h # total number of firms

  # assign exogenous firm-level variables


  #' let X = (PL, PM, PY)
  #' let X_firm = (omega, omegal, K, eta)
}

simAR1 <- function(rho, X) {
  #' simulate an AR1 process forward for a vector of size N_f
  #' @param rho coefficient vector
  #' @param X vector of size N_f to be simulated forward t periods
  #' @param t number of periods to simulate forward



}

simulateOnePeriod <- function(N, N_h) {

  #' @param N number of markets
  #' @param N_h number of firms in market within 2-digit industry
  #'



  # initialize matrix of endogenous variables
  theta0 <- matrix(NA, 12*N_f)
  theta0[1:(8*N_f)] <- rnorm(8*N_f, mean = 0, sd = 1)
  theta0[((8*N_f+1)):(10*N_f)] <- as.matrix(replicate(2*N_f/N_h, matrix(rnorm(n = 1, mean = 0, sd = 1),N_h)))
  theta0[(10*N_f+1):(12*N_f)] <- replicate(2, matrix(rnorm(n = 1, mean = 0, sd = 1), N_f))

  nleqslv(theta0, systemOfEqs, X, X_firm, method = "Broyden", )

}

systemOfEqs <- function(theta, X, X_firm) {
#' Set up system of equations
#'
#' Takes in guess for parameter vector and spits out system of equations to be optimized
#'
#' @param theta guess for parameter vector (in this case vector of exogenous variables)
#' @param beta calibrated parameters
#' @param X exogenous variables
#' let theta = (mu, Pih, lambda, Sih, Lih, Mih, Yih, Yih_star, Yh, Ph, Y, P)
#' let beta = (beta_L, beta_M, v, sigma_CES, sigma, epsilon, N, N_h)
#' let X = (PL, PM, PY)
#' let X_firm = (omega, omegal, K, eta)
#' organize variable matrics as N_f*H (for each t) (so one market, all firms in the market, next market...) (so market-level variables have all subsequent N_h entries equal)



  # assign guesses of endogenous variables (take powers of two to force non-negativity constraints)
  mu <- theta[1:N_f]^2
  Pih <- theta[(N_f+1):(2*N_f)]^2
  lambda <- theta[(2*N_f+1):(3*N_f)]^2
  Sih <- theta[(3*N_f+1):(4*N_f)]^2
  Lih <- theta[(4*N_f+1):(5*N_f)]^2
  Mih <- theta[(5*N_f+1):(6*N_f)]^2
  Yih <- theta[(6*N_f+1):(7*N_f)]^2
  Yih_star <-theta[(7*N_f+1):(8*N_f)]^2
  Yh <- theta[(8*N_f+1):(9*N_f)]^2
  Ph <- theta[(9*N_f+1):(10*N_f)]^2
  Y <-  theta[(10*N_f+1):(11*N_f)]^2
  P <- theta[(11*N_f+1):(12*N_f)]^2



  # assign exogenous variables
  PL <- X[1]
  PM <- X[2]
  PY <- X[3]
  omega <- X_firm[1]
  omegal <- X_firm[2]
  K <- X_firm[3]
  eta <- X_firm[4]


  # set up equations
  r <- matrix(NA, length(theta))

  f <- ((1-beta_L-beta_M)*K^sigma_CES + beta_L*(exp(omegal)*L)^(sigma_CES) + beta_M * M^(sigma_CES))^(v/sigma_CES)
  B <- P_L^(sigma_CES/(sigma_CES-1)) * beta_L^(-1/(sigma_CES-1)) + P_M^(sigma_CES/(sigma_CES-1)) * beta_M^(-1/(sigma_CES-1))

  r[1:N_f] <- Y - sum(Yh[seq(1,N*N_h,N)]) # sum over market outputs
  r[(N_f+1):(2*N_f)] <- Yh - (Ph/P)^(Sigma) * Y
  r[(2*N_f+1):(3*N_f)] <- Yih - (Pih / Ph)^(epsi) * Yh
  r[(3*N_f+1):(4*N_f)] <- Ph - as.matrix(unlist(lapply(unlist(lapply(seq(1,N*N_h,N), function(x) sum(Pih[x:(x+N-1)]^(1-epsi)))),
                      function(x) matrix(x,N_h,1))))
  r[(4*N_f+1):(5*N_f)] <- mu - (epsi/(epsi-1)) / (1 - ((epsi/Sigma)-1)/(epsi-1) *  Sih)
  r[(5*N_f+1):(6*N_f)] <- Sih - (Pih*Yih)(Ph*Yh)
  r[(6*N_f+1):(7*N_f)] <- Yih_star - f*exp(omega)
  r[(7*N_f+1):(8*N_f)] <- PM - lambda*v*f^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*M^(sigma_CES-1)*exp(omega)
  r[(8*N_f+1):(9*N_f)] <- PL - lambda*v*f^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*L^(sigma_CES-1)*exp(omegal)^(sigma_CES)*exp(omega)
  r[(9*N_f+1):(10*N_f)] <- lambda - (1/v)*(beta_L*L^(sigma_CES) + beta_M * M^(sigma_CES))^(1/sigma_CES-1) f^(sigma_CES/v-1)*B^((sigma_CES-1)/sigma_CES)
  r[(10*N_f+1):(11*N_f)] <- P - lambda*mu
  r[(11*N_f+1):(12*N_f)] <- Yih - Yih_star*exp(eta)


}
