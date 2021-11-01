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
  Q = t(Q) %*% Q
  
  phi = param[((m+m*2+m)+(k+2)*(k+2)+1):((m+m*2+m)+(k+2)*(k+2)+(k+2))]
  
  Phi = param[((m+m*2+m)+(k+2)*(k+2)+(k+2)+1):((m+m*2+m)+(k+2)*(k+2)+(k+2)+(k+2)*(k+2))]
  Phi = matrix(Phi, k+2, k+2)
  
  y = as.matrix(t(y[,2:63]))
  names(y) = NULL
  
  # Computes total loglikelihood given a,B,H,Q,phi,Phi and y
  res = KalmanFilter(a, B, H, Q, phi, Phi, y) 
  print((LogLikelihood(res$V, res$u)))
  -(LogLikelihood(res$V, res$u))
}

kalman_optimizer <- function(y) {
  # temp
  m = 60
  k = 2
  
  set.seed(123)
  
  a = runif(m, 0, 0.001)
  B = runif(m*2, 0, 0.001)
  H = runif(m, 0, 1)
  Q = runif((k+2)*(k+2), 0, 1)
  
  phi = rep(0.01, (k+2)) # Lower bound = 0
  Phi = rep(0.01, (k+2)*(k+2)) # Lower bound = 0
  
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
  param = optim(
    par = init_param,
    fn = loss,
    method = "L-BFGS-B",
    lower = rep(0, length(init_param)),
    upper = rep(4, length(init_param)),
    y = y,
    control = list(trace=1)
  )
  
  
  # param = nlminb(init_param, loss, y=y,control=list(trace=1))
  
  # loss(init_param, y)  
  
  # Smoother (?)
  
  param
}
