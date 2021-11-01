
LogLikelihood <- function(V,u) {
  #' Prediciton Error Decomposition representation
  #'
  #' @param V: B P_(t|t-dt) B' + H
  #' @param u: y - (a + B(phi + Phi X_(t-dt)))
  #'
  #' @return loglikelihood
  
  loss = 0
  T = dim(V)[3]
  for (t in 1:T){
    mV = as.matrix(V[,,t])
    loss = loss + (- 0.5 * log(det(mV)) - 0.5 * t(u[,t]) %*% chol2inv(mV) %*% u[,t])
  }
  loss
}
