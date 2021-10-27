
LogLikelihood <- function(V,u) {
  #' Prediciton Error Decomposition representation
  #'
  #' @param V: B P_(t|t-dt) B' + H
  #' @param u: y - (a + B(phi + Phi X_(t-dt)))
  #'
  #' @return loglikelihood
  
  return (- 0.5 * log(det(V)) - 0.5 * t(u) %*% chol2inv(V) %*% u)
  
  # loss = 0
  # for (t in 1:T){
  #   loss = loss + (- 0.5 * log(det(V[,,t])) - 0.5 * t(u[,t]) %*% chol2inv(V[,,t]) %*% u[,t])
  # }
  # return loss
}


# m = 60
# u = rnorm(m+2)
# V = matrix(runif((m+2)^2)*2-1, ncol=(m+2))
# V = t(V) %*% V
# LogLikelihood(V,u)
