loss <- function(param, y) {
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
  
  
  # # Run optimization
  # opt_param = optim(
  #   par = init_param,
  #   fn = loss,
  #   method = "BFGS",
  #   # method = "CG",
  #   y = y,
  #   control = list(trace=1)
  # )

  
  # loss(init_param, y)
  param = nlminb(init_param, loss, y=y, control=list(trace=1))
  # param = nlm(loss,init_param, y=y,steptol=1e-4,gradtol=1e-4, print.level=2)
  # param
}


