initialize_parameters <- function() {
  m = 60
  
  delta_pi = c(0,0,0)
  delta_r = c(0,0,0)
  K = c(0,0,0)
  sigma_pi = c(0,0,0)
  sigma_s = c(0,0,0)
  eta_s = 0
  lambda = c(0,0)
  Lambda = c(0,0,0,0)
  h = rep(0,m)
  
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

define_parameters <- function() {
  library(expm)
  h = 1/12
  
  # (eq. 3.1 Pelsser (2019))
  a = rbind(
    rbind(0,0),
    delta_pi[1] - 0.5 * t(sigma_pi) %*% sigma_pi,
    delta_r[1] + eta_s - 0.5 * t(sigma_s) %*% sigma_s
  )
  
  A = rbind(
    -K, 
    delta_pi[2:3], 
    delta_r[2:3]
  )
  
  C = rbind(
    cbind(diag(2), matrix(0,2,2)),
    as.vector(sigma_pi),
    as.vector(sigma_s)
  )
  
  r = eigen(C)
  U = r$vectors
  D = r$values
  Uinv = solve(U)
  alpha <- function(x) sapply(x, function(x) {if (x == 0) 1 else ((exp(x) - 1) / x)})
  F_ = diag(h * alpha(D * h))
  
  phi = U %*% F_ %*% Uinv %*% a
  Phi = expm(C * h)
  
  # Create V matrix with loop
  V = matrix(0, 4, 4)
  for (i in 1:4) {
    for (j in 1:4) {
      V[i,j] = (Uinv %*% C %*% t(C) %*% t(Uinv))[i,j] * h * alpha((D[i] + D[j]) * h)
    }
  }
  
  Q = U %*% V %*% t(U)
  
  B <- function(tau) {
    sapply(tau, function(tau) if (tau==0) 0 else solve(t(K) + t(Lambda)) %*% expm(-(t(K) + t(Lambda))*tau) %*% delta_r[2:3])
  }
  
  Aprime <- function(tau) {
    sapply(tau, function(tau)  if (tau==0) 0 else -t(B(tau)) %*% lambda + 0.5 * t(B(tau)) %*% B(tau) - delta_r[1])
  }
  
  # Bprime <- function(tau) {
  #   - (t(K) + t(Lambda)) %*% B(tau) - delta_r[2:3]
  # } 
  
  A <- function(tau) {
    sapply(tau, function(tau) integrate(Aprime, 0, tau)$value)
  }
  
  
  
  list(phi=phi, Phi=Phi, Q=Q)
}