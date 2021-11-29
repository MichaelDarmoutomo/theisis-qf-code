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
  Lambda = matrix(lambda, 2, 2)
  h = param[24:length(param)]**2

  l = define_parameters(delta_pi,delta_r,K,sigma_pi,sigma_s,eta_s,lambda, Lambda, h, m)

  # Computes total loglikelihood given a,B,H,Q,phi,Phi and y
  res = KalmanFilter(l$a, l$B, l$H, l$Q, l$phi, l$Phi, y) 
  -(LogLikelihood(res$V, res$u))
}

kalman_optimizer <- function(y) {

  set.seed(123)
  init_param <- initialize_parameters()
  
  lb = rep(-Inf, length(init_param))
  lb[1] = 0
  lb[4] = 0
  lb[24:length(init_param)] = 0

  # Run optimization
  opt_param = optim(
    par = init_param,
    fn = loss,
    method = "L-BFGS-B",
    lower= lb,
    y = y,
    control = list(trace=1)
  )

  # param = nlminb(init_param, loss, y=y, lower=lb,control=list(trace=1))
  # param = nlm(loss,init_param, y=y,steptol=1e-4,gradtol=1e-4, print.level=2)
  # param
}


