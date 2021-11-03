loss <- function(param, processed_param, y) {
  m = 60
  k = 2
  
  delta_pi = param[1:3]
  delta_r = param[4:6]
  K = param[7:9]
  K = matrix(c(K[1],K[2],0,K[3]), 2, 2)
  sigma_pi = c(param[10:12], 0)
  sigma_s = param[13:16]
  eta_s = param[17]
  lambda = param[18:19]
  Lambda = param[20:23]
  Lambda = matrix(lambda, 2, 2)
  h = param[24:length(param)]

  l = define_parameters(delta_pi,delta_r,K,sigma_pi,sigma_s,eta_s,lambda, Lambda, h)
  
  y = as.matrix(t(y[,2:63]))
  names(y) = NULL
  
  # Computes total loglikelihood given a,B,H,Q,phi,Phi and y
  res = KalmanFilter(l$a, l$B, l$H, l$Q, l$phi, l$Phi, y) 
  -(LogLikelihood(res$V, res$u))
}

kalman_optimizer <- function(y) {
  # temp
  m = 60
  k = 2
  
  set.seed(123)
  init_param <- initialize_parameters()
  
  
  # Run optimization
  # opt_param = optim(
  #   par = init_params,
  #   fn = loss,
  #   method = "L-BFGS-B",
  #   # lower = rep(0, length(init_params)),
  #   # upper = rep(1, length(init_params)),
  #   y = y,
  #   control = list(trace=0)
  # )
  # 
  
  param = nlminb(init_param, loss, y=y,lower = rep(0, length(init_param)),control=list(trace=1))

  # param
}


