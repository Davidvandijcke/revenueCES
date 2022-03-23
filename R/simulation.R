#as.matrix(unlist(lapply(test, function(x) matrix(x, 5,1))))

simulation <- function(t = 40, seed = 32424) {

  set.seed(seed)

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

  # assign parameters of exogenous dynamics
  rho_L <- 0.67 # persistence of PL
  sigma_L <- 0.06 # sd of PL
  rho_M <- 0.92 # persistence of PM
  sigma_M <- 0.06 # sd of PM
  rho_D <- 0.86 # persistence of demand
  sigma_D <- 0.68 # sd of demand
  rho_o <- 0.91 # persistence of omega
  sigma_o <- 0.68 # sd of omega
  rho_ol <- 0.92 # persistence of omegal
  sigma_ol <- 0.75 # sd of omegal
  rho_K <- 0.66 # persistence of capital
  sigma_K <- 0.66 # sd of capital
  sigma_eta <- 0.122 # sd measurement error

  # create vectors of exogenous aggregate variables
  X <- exp(rbind(
        repmat(simAR1(rho_L, sigma_L, rnorm(1, mean = 0, sd = sigma_L), t), N_f, 1), # PL
        repmat(simAR1(rho_M, sigma_M, rnorm(1, mean = 0, sd = sigma_M), t), N_f, 1), # PM
        repmat(simAR1(rho_D, sigma_D, rnorm(1, mean = 0, sd = sigma_D), t), N_f, 1) # PM
  )) # take exponent because AR(1) process is for logs

  # create vector of exogenous firm variables
  X_firm <- exp(rbind(
    simAR1(rho_o, sigma_o, rnorm(N_f, mean = 0, sd = sigma_o), t), # omega
    simAR1(rho_ol, sigma_ol, rnorm(N_f, mean = 0, sd = sigma_ol), t), # omegal
    simAR1(rho_K, sigma_K, rnorm(N_f, mean = 0, sd = sigma_K), t), # K
    matrix(rnorm(N_f*t, mean = 0, sd = sigma_eta), nrow = N_f, ncol = t) # eta
  ))




}




simulateOnePeriod <- function(N, N_h, X, X_firm) {

  #' @param N number of markets
  #' @param N_h number of firms in market within 2-digit industry
  #'

  # initialize matrix of endogenous variables
  theta0 <- matrix(NA, 12*N_f) # preassign

  theta0[1:(8*N_f)] <- rnorm(8*N_f, mean = 0, sd = 1) # firm-level variables
  theta0[((8*N_f+1)):(10*N_f)] <- as.matrix(replicate(2*N_f/N, matrix(rnorm(n = 1, mean = 0, sd = 1),N))) # market-level variables
  theta0[(10*N_f+1):(12*N_f)] <- replicate(2, matrix(rnorm(n = 1, mean = 0, sd = 1), N_f)) # aggregate variables


  solved <- nleqslv(theta0, systemOfEqs, jac = NULL, X = X[,1], X_firm[,1], method = "Broyden")

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
  PL <- X[1:(N_f)]
  PM <- X[(N_f+1):(2*N_f)]
  PY <- X[(2*N_f+1):(3*N_f)]
  omega <- X_firm[1:(N_f)]
  omegal <- X_firm[(N_f+1):(2*N_f)]
  K <- X_firm[(2*N_f+1):(3*N_f)]
  eta <- X_firm[(3*N_f+1):(4*N_f)]

  # set up equations
  r <- matrix(NA, length(theta))

  f <- ((1-beta_L-beta_M)*K^sigma_CES + beta_L*(exp(omegal)*Lih)^(sigma_CES) + beta_M * Mih^(sigma_CES))^(v/sigma_CES)
  B <- PL^(sigma_CES/(sigma_CES-1)) * beta_L^(-1/(sigma_CES-1)) + PM^(sigma_CES/(sigma_CES-1)) * beta_M^(-1/(sigma_CES-1))

  r[1:N_f] <- Y - sum(Yh[seq(1,N*N_h,N)]) # CES aggregator over market outputs
  r[(N_f+1):(2*N_f)] <- Yh - (Ph/P)^(Sigma) * Y # CES aggregator over firm outputs
  r[(2*N_f+1):(3*N_f)] <- Yih - (Pih / Ph)^(epsi) * Yh # market output inverse CES share of aggregate
  r[(3*N_f+1):(4*N_f)] <- Ph - as.matrix(unlist(lapply(unlist(lapply(seq(1,N*N_h,N), function(x) sum(Pih[x:(x+N-1)]^(1-epsi)))),
                      function(x) matrix(x,N,1)))) # CES sum over firm prices (need to sum over each N firms in each market first and then expand these again)
  r[(4*N_f+1):(5*N_f)] <- mu - (epsi/(epsi-1)) / (1 - ((epsi/Sigma)-1)/(epsi-1) *  Sih) # markup depends on elasticity and market share
  r[(5*N_f+1):(6*N_f)] <- Sih - (Pih*Yih)/(Ph*Yh) # definition market share
  r[(6*N_f+1):(7*N_f)] <- Yih_star - f*exp(omega) # production function
  r[(7*N_f+1):(8*N_f)] <- PM - lambda*v*f^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Mih^(sigma_CES-1)*exp(omega) # FOC cost min materials
  r[(8*N_f+1):(9*N_f)] <- PL - lambda*v*f^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Lih^(sigma_CES-1)*exp(omegal)^(sigma_CES)*exp(omega) # FOC cost min labor
  r[(9*N_f+1):(10*N_f)] <- lambda - (1/v)*(beta_L*Lih^(sigma_CES) + beta_M * Mih^(sigma_CES))^(1/sigma_CES-1) * f^(sigma_CES/v-1)*B^((sigma_CES-1)/sigma_CES) # marginal cost from cost function
  r[(10*N_f+1):(11*N_f)] <- P - lambda*mu # price equals markup times marginal cost
  r[(11*N_f+1):(12*N_f)] <- Yih - Yih_star*exp(eta) # measurement error on output

  return(r)
}
