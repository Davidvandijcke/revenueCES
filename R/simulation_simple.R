#as.matrix(unlist(lapply(test, function(x) matrix(x, 5,1))))

simulation <- function(t = 40, seed = 32424) {

  set.seed(seed)

  # assign calibrated parameters to global scope
  beta_L <<- 0.4
  beta_M <<- 0.4
  v <<- 0.8 # constant returns to scale
  sigma_CES <<- 0.5 # Raval (2019)
  Sigma <<- 1.1 # as in De Ridder et. al.
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
    repmat(simAR1(rho_L, sigma_L, rnorm(1, mean = log(100), sd = sigma_L), t, b0 = log(100)), N_f, 1), # PL
    repmat(simAR1(rho_M, sigma_M, rnorm(1, mean = log(100), sd = sigma_M), t, b0 = log(100)), N_f, 1), # PM
    repmat(simAR1(rho_D, sigma_D, rnorm(1, mean = log(7), sd = sigma_D), t, b0 = log(7)), N_f, 1) # PY
  )) # take exponent because AR(1) process is for logs

  # create vector of exogenous firm variables
  X_firm <- exp(rbind(
    simAR1(rho_o, sigma_o, rnorm(N_f, mean = log(1.2), sd = sigma_o), t, b0 = log(1.2)), # omega
    simAR1(rho_ol, sigma_ol, rnorm(N_f, mean = log(1.2), sd = sigma_ol), t, b0 = log(1.2)), # omegal
    simAR1(rho_K, sigma_K, rnorm(N_f, mean = log(3), sd = sigma_K), t, b0 =  log(3)), # K
    matrix(rnorm(N_f*t, mean = log(1), sd = sigma_eta), nrow = N_f, ncol = t) # eta
  ))




}




simulateOnePeriod <- function(N, N_h, X, X_firm) {

  #' @param N number of markets
  #' @param N_h number of firms in market within 2-digit industry
  #'

  # initialize matrix of endogenous variables
  theta0 <- matrix(NA, 12*N_f) # preassign

  theta0[((N_f+1)):(2*N_f)] <- rnorm(N_f, mean = log(100), sd = 0.06) # firm-level log prices
  theta0[((4*N_f+1)):(6*N_f)] <-  rnorm(N_f, mean = log(100), sd = 0.1) # firm-level log inputs
  theta0[((9*N_f+1)):(10*N_f)] <- as.matrix(replicate(N_f/N, matrix(rnorm(n = 1, mean = log(100), sd = 0.06),N))) # market-level log prices
  theta0[(11*N_f+1):(12*N_f)] <- replicate(1, matrix(rnorm(n = 1, mean = log(100), sd = 0.05), N_f)) # aggregate log prices
  theta0[(10*N_f+1):(11*N_f)] <- log(X_firm[(2*N_f+1):(3*N_f),1]) - Sigma*(theta0[(11*N_f+1):(12*N_f)]) # aggregate output = aggregate revenue / prices^Sigma
  theta0[((8*N_f+1)):(9*N_f)] <- -Sigma * ( theta0[((9*N_f+1)):(10*N_f)] - theta0[(11*N_f+1):(12*N_f)] ) +  theta0[(10*N_f+1):(11*N_f)]  # market-level output
  theta0[((6*N_f+1)):(7*N_f)] <- -epsi*(theta0[((N_f+1)):(2*N_f)] - theta0[(11*N_f+1):(12*N_f)]) + theta0[((8*N_f+1)):(9*N_f)] # firm-level output
  theta0[((7*N_f+1)):(8*N_f)] <- theta0[((6*N_f+1)):(7*N_f)] + log(X_firm[(3*N_f+1):(4*N_f), 1])  # firm-level output with measurement error
  theta0[((3*N_f+1)):(4*N_f)] <- rnorm(N_f, mean = log(0.005), sd = 0.1) # firm market share
  theta0[1:(N_f)] <-  log((epsi/(epsi-1))) - log((1 - ((epsi/Sigma)-1)/(epsi-1) *  exp(theta0[((3*N_f+1)):(4*N_f)]) ) ) # firm-level markups
  theta0[((2*N_f+1)):(3*N_f)] <- theta0[1:(N_f)] - theta0[((N_f+1)):(2*N_f)]   # firm-level log marginal costs = log markup - log price



  theta <- repmat(log(1), N_f)

  #rnorm(N_f, mean = log(0.005), sd = 0.1) # firm market share = firm revenue / market revenue


  solved <- nleqslv(theta0, systemOfEqs, jac = NULL, X = X[,1], X_firm[,1], method = "Broyden")

  solved <- optim(as.vector(theta0), systemOfEqs, X = X[,1], X_firm = X_firm[,1])
  solved <- BB::BBsolve(as.vector(theta0), systemOfEqs, X = X[,1], X_firm = X_firm[,1])

  solved <- pracma::fsolve(x = as.vector(theta0), f = systemOfEqs, X = X[,1], X_firm = X_firm[,1])


}
systemOfEqs <- function(gamma, X, X_firm) {
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
  theta <- exp(gamma) # impose non-negativity/
  Pih <- theta[(1:(N_f))] # Ph
  Lih <- theta[((N_f+1)):(2*N_f)]
  Mih <- theta[((2*N_f+1)):(3*N_f)]
  # assign exogenous variables
  PL <- X[1:(N_f)]
  PM <- X[(N_f+1):(2*N_f)]
  PY <- X[(2*N_f+1):(3*N_f)]
  omega <- X_firm[1:(N_f)]
  omegal <- X_firm[(N_f+1):(2*N_f)]
  K <- X_firm[(2*N_f+1):(3*N_f)]
  eta <- X_firm[(3*N_f+1):(4*N_f)]
  # intermediate equations
  Ph <- as.matrix(unlist(lapply(unlist(lapply(seq(1,N*N_h,N), function(x) sum(Pih[x:(x+N-1)]^(1-epsi))^(1/(1-epsi)) )),
                                function(x) matrix(x,N,1)))) # CES sum over firm prices (need to sum over each N firms in each market first and then expand these again)
  Yh <- (Ph)^(-Sigma) * PY
  Yih <- (Pih / Ph)^(-epsi) * Yh # market output inverse CES share of aggregate
  Y <- sum(Yh[seq(1,N*N_h,N)]^((Sigma-1)/Sigma))^(Sigma/(Sigma-1))
  P <- (Yh/Y)^(1/Sigma)*Ph
  #P <- (PY / Y)^(1/Sigma)
  Sih <- (Pih*Yih)/(Ph*Yh) # definition market share
  mu <- (epsi/(epsi-1)) / (1 - ((epsi/Sigma)-1)/(epsi-1) *  Sih) # markup
  lambda <- Pih / mu # marginal cost
  f <- ((1-beta_L-beta_M)*K^(sigma_CES) + beta_L*(omegal*Lih)^(sigma_CES) + beta_M * Mih^(sigma_CES))^(v/sigma_CES) # production technology
  B <- PL^(sigma_CES/(sigma_CES-1)) * beta_L^(-1/(sigma_CES-1)) + PM^(sigma_CES/(sigma_CES-1)) * beta_M^(-1/(sigma_CES-1))
  # set up equations
  r <- matrix(NA, length(theta))
  r[(1):(N_f)]  <- PL - lambda*v*(f)^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Lih^(sigma_CES-1)*omegal^(sigma_CES)*omega # FOC cost min labor
  r[(N_f+1):(2*N_f)] <- PM - lambda*v*(f)^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Mih^(sigma_CES-1)*omega # FOC cost min materials
  r[(2*N_f+1):(3*N_f)] <- PY - P^(Sigma)*Y
  #Yih -  f*omega
  print(summary(r))
  #return(error^2)
  return(as.vector(r))
}
