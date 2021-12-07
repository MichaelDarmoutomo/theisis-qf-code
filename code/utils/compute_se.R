compute_se <- function(param, data) {
  Hess = maxLik::numericHessian(function(p) loss(p, data), t0=param)
  negInvHess = solve(Hess)
  sqrt(diag(negInvHess))
}