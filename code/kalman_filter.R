loss <- function(a, B, H, Q, phi, Phi, y) {
  # Computes total loglikelihood given a,B,H,Q,phi,Phi and y
  res = kalman(a, B, H, Q, phi, Phi, y) 
  return -(loglikelihood(res$V, res$u))
}

kalman_optimizer <- function(y) {
  
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
  param = optim(init_param, loss)
  
  
  # Smoother (?)
  
  return (param)
}
