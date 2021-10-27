KalmanFilter <- function(parameter_vector, y)
# This function runs the Kalman filter for the scalar AR(1) model plus
# noise with a diffuse prior (roughly)

# % Extract lenght of the data
T = dim(y)[-1]

# % Extract the stuff we need from the input arguments
F = parameter_vector(1,1)
Q = parameter_vector(2,1)
R = parameter_vector(3,1)

# % The Kalman filter for AR(1)+noise
for (t in 1:T) {
  
  # % Prediction for t=1
  if (t == 1) {
    Xhat[t] = 0;
    predictedP[t] = 10^6;
  } 
  else {
    # % Prediction step for any other t
    Xhat[t] = phi + Phi * X(t-1) # Predicted state
    Phat[t] = F * P(t-1) * t(F) + Q # Predicted covariance (Q is variance of state eq.)
  }
  # % Update for any t (R is variance of meas. eq.)
  # X[t] = predictedX[t] + predictedP[t] * 1/( predictedP[t] + R ) * ( y[t] - predictedX[t] )
  # P[t] = predictedP[t] - predictedP[t] * 1/( predictedP[t] + R ) * predictedP[t]
  K[t] = Phat %*% t(B) %*% chol2inv(Phat[t])
  # State update
  X[t] = Xhat[t] + K[t] %*% (y[t] - (a + B*Xhat[t]))

  
  }


