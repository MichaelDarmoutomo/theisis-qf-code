initialize_parameters <- function() {
  m = 6

  # delta_pi = rnorm(3, 0 , 0.1) # c(0,0,0)
  # delta_r = rnorm(3, 0 , 0.1) # c(0,0,0)
  # K = runif(3, 0, 1) # c(0,0,0)
  # sigma_pi = rnorm(3, 0, 0.2) # c(0,0,0)
  # sigma_s = rnorm(4, 0, 0.2) #c(0,0,0,0)
  # eta_s = rnorm(1, 0, 0.2)
  # lambda = rnorm(2, 0, 0.1) # c(0,0)
  # Lambda = rnorm(4, 0, 0.1) # c(0,0,0,0)
  # h = rep(0.001,m)
  
  # delta_pi = c(0.02,-0.001, -0.001)
  # delta_r =  c(0.01, -0.001, -0.002)
  # K = c(0.5, 0.5, 0.2)
  # sigma_pi = c(0.1,0.1,0.1)
  # sigma_s = c(0.1,0.1,0.1,0.1)
  # eta_s = 0.05
  # lambda =  c(0.5,0.1)
  # Lambda =  c(0.2,0.4,-0.5,-1)
  # h = rep(0.01,m)
  
  delta_pi = c(0.0158, -0.0028, -0.0014)
  delta_r = c(0.0097, -0.0094, -0.0024)
  K = c(0.0479, 0.5440, 1.2085)
  sigma_pi = c(-0.0010, 0.0013, 0.0055)
  sigma_s = c(-0.0483, 0.0078, 0.0010, 0.1335)
  eta_s = 0.0451
  lambda = c(0.6420, -0.0240)
  Lambda = c(0.1710, 0.3980, -0.5140, -1.1470)
  h = c(0.0038, 0.0003, 0.0003, 0.0000, 0.0008, 0.0021)

  c(delta_pi,
    delta_r,
    K,
    sigma_pi,
    sigma_s,
    eta_s,
    lambda, 
    Lambda,
    h)
}

define_parameters <- function(delta_pi,delta_r,K,sigma_pi,sigma_s,eta_s,lambda, Lambda, h_, m) {
  h = 1/12
  
  # (eq. 3.1 Pelsser (2019))
  a_ = rbind(
    rbind(0,0),
    delta_pi[1] - 0.5 * t(sigma_pi) %*% sigma_pi,
    delta_r[1] + eta_s - 0.5 * t(sigma_s) %*% sigma_s
  )
  
  A_ = cbind(
    rbind(
    -K, 
    delta_pi[2:3], 
    delta_r[2:3]
    ),
    matrix(0, nrow=4, ncol=2)
  )
  
  C_ = rbind(
    cbind(diag(2), matrix(0,2,2)),
    as.vector(sigma_pi),
    as.vector(sigma_s)
  )
  
  r = eigen(A_)
  U = r$vectors
  D = r$values
  Uinv = solve(U)
  alpha <- function(x) ifelse(x == 0, 1, ((exp(x) - 1) / x))
  F_ = diag(h * alpha(D * h))
  
  phi = U %*% F_ %*% Uinv %*% a_
  Phi = expm(A_ * h)
  
  # Create V matrix with loop
  V = matrix(0, 4, 4)
  for (i in 1:4) {
    for (j in 1:4) {
      V[i,j] = (Uinv %*% C_ %*% t(C_) %*% t(Uinv))[i,j] * h * alpha((D[i] + D[j]) * h)
    }
  }
  
  Q = U %*% V %*% t(U)
  
  fB <- function(tau) {
    M = t(K + Lambda)
    Minv = solve(M)
    sapply(tau, function(tau) if (tau==0) 0 else Minv %*% (expm(-tau * M) - diag(2)) %*% delta_r[2:3])
  }
  
  fAprime <- function(tau) {
    sapply(tau, function(tau)  if (tau==0) 0 else -t(fB(tau)) %*% lambda + 0.5 * t(fB(tau)) %*% fB(tau) - delta_r[1])
  }
  
  # Bprime <- function(tau) {
  #   - (t(K) + t(Lambda)) %*% B(tau) - delta_r[2:3]
  # } 
  
  fA <- function(tau) {
    # sapply(tau, function(tau) integrate(fAprime, 0, tau)$value)
    tau = c(0, tau)
    p = sapply(seq(length(tau)-1), function(k) integrate(fAprime, tau[k], tau[k+1], abs.tol=0.01)$value)
    cumsum(p)
  }
  
  a_temp = fA(c(1,5,10,15,20,30)) / c(1,5,10,15,20,30)
  a = c(-a_temp, 0, 0)
  
  B = matrix(0, m+2, 4)
  B_temp = t(fB(c(1,5,10,15,20,30))) / c(1,5,10,15,20,30)
  B[1:m,1:2] = -B_temp
  B[(m+1):(m+2), 3:4] = diag(2)

  H = matrix(0, m+2, m+2)
  H[1:m, 1:m] = diag(h_)
  list(a=a, B=B, H=H, Q=Q, phi=phi, Phi=Phi)
}