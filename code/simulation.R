simulate_economy <- function(paramfile, dt, T, nSim, w_0, maturities, savepath=FALSE, parallel=FALSE, save_economy=FALSE, save_path="") {
  library(abind)
  library(expm)
  
  # load('results/Parameters.Rdata')
  load(paramfile)
  # delta_pi = c(0.0158, -0.0028, -0.0014)
  # delta_r = c(0.0097, -0.0094, -0.0024)
  # K = c(0.0479, 0.5440, 0.7085)
  # sigma_pi = c(-0.0010, 0.0013, 0.0055)
  # sigma_s = c(-0.0483, 0.0078, 0.0010, 0.1335)
  # eta_s = 0.0451
  # lambda = c(0.6420, -0.0240)
  # Lambda = c(0.1710, 0.3980, -0.5140, -1.1470)
  # h = c(0.0038, 0.0003, 0.0003, 0.0000, 0.0008, 0.0021)
  # 
  # param=c(delta_pi,
  #   delta_r,
  #   K,
  #   sigma_pi,
  #   sigma_s,
  #   eta_s,
  #   lambda, 
  #   Lambda,
  #   h)
  
  # Define parameters
  delta_pi = param[1:3]
  delta_r = param[4:6]
  K = param[7:9]
  K = matrix(c(K[1],K[2],0,K[3]), 2, 2)
  sigma_pi = c(param[10:12], 0)
  sigma_s = param[13:16]
  eta_s = param[17]
  lambda_ = param[18:19]
  lambda = c(lambda_, 0, t(-sigma_s[1:2]) %*% lambda_[1:2] / sigma_s[4])
  Lambda_ = param[20:23]
  Lambda_ = matrix(Lambda_, 2, 2, byrow=TRUE)
  Lambda = rbind(
    Lambda_, 
    c(0,0))
  Lambda = rbind(Lambda,
                 c(t(-sigma_s[1:3]) %*% Lambda[1:3,1] / sigma_s[4],
                   t(-sigma_s[1:3]) %*% Lambda[1:3,2] / sigma_s[4]))
  h = param[24:length(param)]
  sigma_x = rbind(diag(2), matrix(0, 2, 2))
  
  #Create inflation, stock market, yields
  dZ          = array(rnorm(4*T*nSim, mean=0, sd=dt), c(4, nSim, T))
  Lambda_t    = array(0, c(4, nSim, T))
  r           = array(0, c(T,nSim))
  pi          = array(0, c(T,nSim))
  dPi_        = array(0, c(T,nSim))
  Pi          = array(0, c(nSim,T+1))
  dS_         = array(0, c(nSim,T))
  S           = array(0, c(nSim,T+1))
  dX          = array(0, c(2, nSim, T))
  X           = array(0, c(2, nSim, T+1))
  
  X[,,1] = matrix(1, 2, nSim) * c(2.364, 1.005)
  Pi[,1] = rep(1, nSim)
  S[,1] = rep(1, nSim)
  
  
  for (t in 1:T) {
    r[t,] = delta_r[1] + t(delta_r[2:3]) %*% X[,,t]
    pi[t,] = delta_pi[1] + t(delta_pi[2:3]) %*% X[,,t]
    Lambda_t = lambda + Lambda %*% X[,,t]
    dPi_[t,] = pi[t,] * dt + t(sigma_pi) %*% dZ[,,t]
    Pi[,t+1] =  Pi[,t] + Pi[,t] * t(dPi_[t,])
    dS_[,t] = (r[t,] + eta_s) * dt + t(sigma_s) %*% dZ[,,t]
    S[,t+1] = S[,t] + S[,t] * dS_[,t]
    dX[,,t] = -K %*% X[,,t] * dt + t(sigma_x)  %*% dZ[,,t]
    X[,,t+1] = X[,,t] + dX[,,t]
  }
  
  Pi   = t(Pi[,-dim(Pi)[2]])
  dS_  = t(dS_)
  S    = t(S[,-dim(S)[2]])
  X    = X[,,-dim(X)[3]]
  
  # Bond prices
  fB <- function(tau) {
    sapply(tau, function(tau) if (tau==0) c(0,0) else solve(t(K + Lambda_)) %*% (expm(-tau * (t(K) + t(Lambda_))) - diag(2)) %*% delta_r[2:3])
  }
  
  fAprime <- function(tau) {
    sapply(tau, function(tau)  if (tau==0) 0 else -t(fB(tau)) %*% lambda_ + 0.5 * t(fB(tau)) %*% fB(tau) - delta_r[1])
  }
  
  
  fA <- function(tau) {
    # if tau is consecutive numbers
    # if (sum((tau - shift(tau, 1))[2:length(tau)] == 1) == length(tau) - 1) {
    #   p = sapply(seq(length(tau)-1), function(k) integrate(fAprime, tau[k], tau[k+1])$value)
    #   cumsum(p)
    # } else {
    #   sapply(tau, function(tau) integrate(fAprime, 0, tau)$value)
    # }
    tau = c(0, tau)
    p = sapply(seq(length(tau)-1), function(k) integrate(fAprime, tau[k], tau[k+1])$value)
    cumsum(p)
  }
  
  
  y = array(0, c(length(maturities), nSim, T+1))
  
  if (parallel) {
    library(parallel)
    # maturities_ = maturities[maturities==0] = 1
    
    numCores <- detectCores()
    clust <- makeCluster(numCores-6)
    print(maturities)
    yield_wrapper = function(t) {
      # exp((fA(maturities) / maturities_) + (t(fB(maturities)) / maturities_) %*% X[,,t])
      exp((fA(maturities)) + (t(fB(maturities))) %*% X[,,t])
    }
    clusterExport(clust, varlist=c("X", "expm", "fA"), envir=environment())
    
    on.exit(stopCluster(clust))
    
    print(paste("Starting at",Sys.time()))
    y = parLapply(clust, 1:T, yield_wrapper)
    print(paste("Stopping at",Sys.time()))
    
  } else {
    print(paste("Starting at",Sys.time()))
    for (t in 1:T) {
      # maturities_ = maturities[maturities==0] = 1
      y[,,t] = exp((fA(maturities)) + (t(fB(maturities))) %*% X[,,t])
    }
    print(paste("Stopping at",Sys.time()))
  }

  # data = abind(y,
  w = w_0 * Pi

  P_tmp = array(unlist(y), c(length(maturities), nSim, T+1))
  P1 = P_tmp[1:65,,]
  P2 = P_tmp[2:66,,]
  
  e = list(
    delta_r = delta_r,
    delta_pi = delta_pi,
    dZ = dZ,
    Lambda = Lambda_t,
    r = r,
    pi = pi,
    dPi_ = dPi_,
    Pi = Pi,
    dS_ = dS_,
    S = S,
    dX = dX,
    X = X,
    P = P1,
    P2 = P2,
    w = w
  )            
  # array(log(Pi), c(1,100,214)),
  #              array(log(S), c(1,100,214)), along=1)
  # 

  if (save_economy) {
    save(e, file=save_path)
  }
  return(e)
}







