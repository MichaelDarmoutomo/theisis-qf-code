loss <- function(param, y) {
  m = dim(y)[1]-2
  
  delta_pi = param[1:3]
  delta_r = param[4:6]
  K = param[7:9]
  K = matrix(c(K[1],K[2],0,K[3]), 2, 2)
  sigma_pi = c(param[10:12], 0)
  sigma_s = param[13:16]
  eta_s = param[17]
  lambda = param[18:19]
  Lambda = param[20:23]
  Lambda = matrix(Lambda, 2, 2, byrow=TRUE)
  h = param[24:length(param)]**2

  l = define_parameters(delta_pi,delta_r,K,sigma_pi,sigma_s,eta_s,lambda, Lambda, h, m)

  res = KalmanFilter(l$a, l$B, l$H, l$Q, l$phi, l$Phi, y) 
  
  -(LogLikelihood(res$V, res$u))
}

kalman_optimizer <- function(y, maxit=10000) {

  init_param <- initialize_parameters()
   
  # opt_param = optim(
  #   par = init_param,
  #   fn = loss,
  #   y = y,
  #   hessian = TRUE,
  #   control = list(trace=2, maxit=maxit)
  # )
  
  opt_param = trust.optim(
    init_param,
    function(p) loss(p, y),
    function(p_) maxLik::numericGradient(function(p) loss(p, data), p_),
    # function(p_) numDeriv::grad(function(p) loss(p, data), p_),
    control=list(maxit=maxit)
  )
  
  
  opt_param
}


