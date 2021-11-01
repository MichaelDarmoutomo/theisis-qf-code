loss <- function(param, y) {
  m = 60
  k = 2
  
  # Prepare variables a, B, H, Q, phi, Phi (with restrictions)
  a = rep(0, m+2)
  a[1:m] = param[1:m]
  
  B = matrix(0, m+2, k+2)
  B[(1:m),(1:k)] = matrix(param[(m+1):(m + m*2)], m, k)
  B[(m+1):(m+2), 3:4] = diag(2)
  
  H = matrix(0, m+2, m+2)
  H[1:m, 1:m] = diag(param[(m+m*2+1):(m+m*2+m)])
  
  Q = param[(m+m*2+m+1):((m+m*2+m)+(k+2)*(k+2))]
  Q = matrix(Q, k+2, k+2)
  
  phi = param[((m+m*2+m)+(k+2)*(k+2)+1):((m+m*2+m)+(k+2)*(k+2)+(k+2))]
  
  Phi = param[((m+m*2+m)+(k+2)*(k+2)+(k+2)+1):((m+m*2+m)+(k+2)*(k+2)+(k+2)+(k+2)*(k+2))]
  Phi = matrix(Phi, k+2, k+2)
  
  y = as.matrix(t(y[,2:63]))
  names(y) = NULL
  
  # Computes total loglikelihood given a,B,H,Q,phi,Phi and y
  res = KalmanFilter(a, B, H, Q, phi, Phi, y) 
  -(LogLikelihood(res$V, res$u))
}

kalman_optimizer <- function(y) {
  # temp
  m = 60
  k = 2
  
  a = rnorm(m)/100
  B = rnorm(m*2)/100
  H = runif(m, 0, 1)
  Qtmp = matrix(rnorm((k+2)*(k+2)), k+2, k+2)
  Q = as.vector(t(Qtmp) %*% Qtmp)
  
  phi = rep(1, (k+2))/100
  Phi = rep(1, (k+2)*(k+2))/100
  
  # Initialize (m is # maturities)
  init_param = c(
    a, # (m + 2)
    B, # (m + 2) x (k + 2)
    H, # (m + 2) x (m + 2)
    Q, # (k + 2) x (k + 2)
    phi, # (k + 2) x 1
    Phi  # (k + 2) x (k + 2)
  )
  
  # Run optimization (use some kind of optimizer)
  # [a, B, H, Q, phi, Phi], loss = optimize(reward)
  # param = optim(par=init_param, method="L-BFGS-B", fn=loss, y=y, control=list(trace=1))
  param = nlminb(init_param, loss, y=y,control=list(trace=1))
  
  
  # Smoother (?)
  
  return (param)
}
