
simulation <- function(t = 40, seed = 2335) {

  set.seed(seed)

  beta_L <<- 0.3
  beta_M <<- 0.3
  v <<- 1 # constant returns to scale
  sigma_CES <<- 0.5 # Raval (2019)
  Sigma <<- 1.1 # as in De Ridder et. al.
  epsi <<- 5 # as in De Ridder et. al.
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

  # K <- system.file("extdata", "orbis_sim_sample.csv", packages = "revenueCES", mustWork = TRUE)
  K_mat <- read.csv('/home/antonvocalis/Dropbox (University of Michigan)/revenueCES/inst/extdata/orbis_sim_sample.csv')
  K_mat <- K_mat / 10^10

  D_init <<- sum(K_mat$R)
  P_init <- 1

  # create vectors of exogenous aggregate variables
  X <- exp(rbind(
    repmat(simAR1(rho_L, sigma_L, rnorm(1, mean = log(P_init), sd = sigma_L), t), N_f, 1), # PL
    repmat(simAR1(rho_M, sigma_M, rnorm(1, mean = log(P_init), sd = sigma_M), t), N_f, 1), # PM
    repmat(simAR1(rho_D, sigma_D, rnorm(1, mean = log(D_init), sd = sigma_D), t, b0 = log(D_init)), N_f, 1) # PY
  )) # take exponent because AR(1) process is for logs
  # create vector of exogenous firm variables
  # create vector of exogenous firm variables
  X_firm <- exp(rbind(
    simAR1(rho_o, sigma_o, rnorm(N_f, mean = log(1), sd = sigma_o), t, b0 = log(1)), # omega
    simAR1(rho_ol, sigma_ol, rnorm(N_f, mean = log(1), sd = sigma_ol), t, b0 = log(1)), # omegal
    #simAR1(rho_K, sigma_K, rnorm(N_f, mean = log(K_init), sd = sigma_K), t, b0 = log(K_init)) # K    simAR1(rho_K, sigma_K, rnorm(N_f, mean = log(K_init), sd = sigma_K), t, b0 = log(K_init)) # K
    simAR1(rho_K, sigma_K, K_mat$K, t) # K    simAR1(rho_K, sigma_K, rnorm(N_f, mean = log(K_init), sd = sigma_K), t, b0 = log(K_init)) # K
  ))

  # initialize matrix of endogenous variables
  theta0 <- matrix(NA, 2*N_f) # preassign
  theta0[(1):(2*N_f)] <-  rnorm(N_f, mean = log(1), sd = sigma_K)  # firm level inputs


  for (step in 1:t) { # for each time period
    print(paste0("PERIOD", step))
    solved <- BB::BBsolve(as.vector(theta0), systemOfEqs, X = X[,step], X_firm = X_firm[,step],
                          control = list(maxit = 100000, noimp = 1000, tol = 1.e-7)) # solve system of non-linear equations
    #solved <- pracma::fsolve(x = as.vector(theta0), f = systemOfEqs, X = X[,step], X_firm = X_firm[,step])

    print(solved$message)
    add <- 0.1
    while (solved$convergence > 0) {
      theta0[(1):(2*N_f)] <-  rnorm(N_f, mean = log(1*add), sd = sigma_K)  # firm level inputs
      print("Try again")
      solved <- BB::BBsolve(as.vector(theta0), systemOfEqs, X = X[,step], X_firm = X_firm[,step],
                            control = list(maxit = 2000, noimp = 300, tol = 1.e-7))
      add <- add*10
    }
    df_temp <-  processSimData(solved$par, X = X[,step], X_firm = X_firm[,step],
                               t = step)
    if (step == 1) { df <- df_temp
    } else { df <- rbind(df, df_temp) }

    theta0 <- solved$par # use converged vector from last period as new initial vector
  }
  return(df)


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
  Lih <- theta[(1:(N_f))] # Ph
  Mih <- theta[((N_f+1)):(2*N_f)]

  # assign exogenous variables
  PL <- X[1:(N_f)]
  PM <- X[(N_f+1):(2*N_f)]
  PY <- X[(2*N_f+1):(3*N_f)]
  omega <- X_firm[1:(N_f)]
  omegal <- X_firm[(N_f+1):(2*N_f)]
  K <- X_firm[(2*N_f+1):(3*N_f)]

  # intermediate equations
  f <- ((1-beta_L-beta_M)*K^(sigma_CES) + beta_L*(omegal*Lih)^(sigma_CES) + beta_M * Mih^(sigma_CES))^(v/sigma_CES) # production technology
  Yih <- f*omega
  Yh <- as.matrix(unlist(lapply(unlist(lapply(seq(1,N*N_h,N), function(x) sum(Yih[x:(x+N-1)]^((epsi-1)/epsi))^(epsi/(epsi-1)) )),
                               function(x) matrix(x,N,1))))
  Y <- sum(Yh[seq(1,N*N_h,N)]^((Sigma-1)/Sigma))^(Sigma/(Sigma-1))
  Ph <- (Yh / PY)^(-1/Sigma)
  Pih <- (Yih / Yh)^(-1/epsi) * Ph
  Sih <- (Pih*Yih)/(Ph*Yh) # definition market share
  mu <- (epsi/(epsi-1)) / (1 - ((epsi/Sigma)-1)/(epsi-1) *  Sih) # markup
  lambda <- Pih / mu # marginal cost

  # set up equations
  r <- matrix(NA, length(theta))
  r[(1):(N_f)]  <- PL - lambda*v*(Yih/omega)^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Lih^(sigma_CES-1)*omegal^(sigma_CES)*omega # FOC cost min labor
  r[(N_f+1):(2*N_f)] <- PM - lambda*v*(Yih/omega)^((v-sigma_CES)/sigma_CES)*(1-beta_L-beta_M)*Mih^(sigma_CES-1)*omega # FOC cost min materials

  #r[(3*N_f+1):(4*N_f)] <- Yih - f*omega

  #Yih -  f*omega
  #print(summary(r))
  #return(error^2)
  return(as.vector(r))
}


processSimData <- function(gamma, X, X_firm, t) {

  theta <- exp(gamma) # impose non-negativity/
  Lih <- theta[(1:(N_f))] # Ph
  Mih <- theta[((N_f+1)):(2*N_f)]
  # assign exogenous variables
  PL <- X[1:(N_f)]
  PM <- X[(N_f+1):(2*N_f)]
  PY <- X[(2*N_f+1):(3*N_f)]
  omega <- X_firm[1:(N_f)]
  omegal <- X_firm[(N_f+1):(2*N_f)]
  K <- X_firm[(2*N_f+1):(3*N_f)]


  # intermediate equations
  f <- ((1-beta_L-beta_M)*K^(sigma_CES) + beta_L*(omegal*Lih)^(sigma_CES) + beta_M * Mih^(sigma_CES))^(v/sigma_CES) # production technology
  Yih <- f*omega
  Yh <- as.matrix(unlist(lapply(unlist(lapply(seq(1,N*N_h,N), function(x) sum(Yih[x:(x+N-1)]^((epsi-1)/epsi))^(epsi/(epsi-1)) )),
                                function(x) matrix(x,N,1))))
  Y <- sum(Yh[seq(1,N*N_h,N)]^((Sigma-1)/Sigma))^(Sigma/(Sigma-1))
  Ph <- (Yh / PY)^(-1/Sigma)
  Pih <- (Yih / Yh)^(-1/epsi) * Ph
  Sih <- (Pih*Yih)/(Ph*Yh) # definition market share
  mu <- (epsi/(epsi-1)) / (1 - ((epsi/Sigma)-1)/(epsi-1) *  Sih) # markup
  lambda <- Pih / mu # marginal cost


  df <- data.frame(period = t, Pih = Pih, Lih = Lih, Mih = Mih, PL = PL, PM = PM,
             PY = PY, omega = omega, omegal = omegal, K = K, Ph = Ph, Yh = Yh,
             Yih = Yih, Y = Y, mu = mu, Sih = Sih)
  return(df)

}
