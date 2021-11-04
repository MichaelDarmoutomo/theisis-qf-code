KalmanFilter <- function(a, B, H, Q, phi, Phi, y) {
  # Set prior
  T = dim(y)[-1]
  X0 = matrix(0, 4, 1)
  P0 = diag(4)
  
  # Define empty arrays and matrices
  Xhat = X = matrix(0, 2+2, T)
  Phat = P = array(0, c(4, 4, T))
  V = array(0, c(dim(H),T))
  u = array(0, dim(y))
  K = array(0, c(dim(t(B)),T))
  
  # Kalman filter
  for (t in 1:T) {
    
    # % Prediction for t == 1
    if (t == 1) {
      Xhat[,t] = X0
      Phat[,,t] = P0
    } 
    else {
      # % Prediction step for any other t
      Xhat[,t] = phi + Phi %*% X[,t-1] # Predicted state
      Phat[,,t] = Phi %*% P[,,(t-1)] %*%  t(Phi) + Q # Predicted covariance (Q is variance of state eq.)
    }
    
    # Variables needed
    V[,,t] = B %*% Phat[,,t] %*% t(B) + H
    u[,t] = y[,t] - (a + B %*% Xhat[,t])
    K[,,t] = Phat[,,t] %*% t(B) %*% solve(V[,,t])
    
    # State update
    X[,t] = Xhat[,t] + K[,,t] %*% u[,t]
    
    # P update
    P[,,t] = (diag(4) - (K[,,t] %*% B)) %*% Phat[,,t]
  }
  return (list(V=V, u=u))
}