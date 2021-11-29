# load parameters
rm(list=ls())
load('results/Parameters.Rdata')

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
Lambda_ = matrix(Lambda_, 2, 2)
Lambda = rbind(
  Lambda_, 
  c(0,0))
Lambda = rbind(Lambda,
  c(t(-sigma_s[1:3]) %*% Lambda[1:3,1] / sigma_s[4],
    t(-sigma_s[1:3]) %*% Lambda[1:3,2] / sigma_s[4]))
h = param[24:length(param)]
sigma_x = rbind(diag(2), matrix(0, 2, 2))

#Create inflation, stock market, yields
dt = 1/12
T = 213
nSim = 100

dZ          = array(rnorm(4*T*nSim, 0, dt), c(4, nSim, T))
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

# Bond prices
fB <- function(tau) {
  sapply(tau, function(tau) if (tau==0) 0 else chol2inv(t(K + Lambda_)) %*% (expm(-tau * (t(K) + t(Lambda_))) - diag(2)) %*% delta_r[2:3])
}

fAprime <- function(tau) {
  sapply(tau, function(tau)  if (tau==0) 0 else -t(fB(tau)) %*% lambda_ + 0.5 * t(fB(tau)) %*% fB(tau) - delta_r[1])
}

Bprime <- function(tau) {
  - (t(K) + t(Lambda_)) %*% B(tau) - delta_r[2:3]
}

fA <- function(tau) {
  p = sapply(seq(length(tau)-1), function(k) integrate(fAprime, tau[k], tau[k+1], rel.tol=.Machine$double.eps^0.18)$value)
  cumsum(p)
}

# exp(fA(0:6) + fB(1) * X[,,,])
