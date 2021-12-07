adjustment_factor <- function(option=1, Tw, Tp) {
  if (option == 1) {
    
    af = function() rep(1, Tw + Tp)
  
  } else if (option == 2) {
    
    af = function() c(rep(3, 30), rep(2, 10), rep(1, 25)
    )
    
  } else if (option == 2) {
    
    af = function(beta) exp(beta[1] * ((Tw + Tp - 1):0) + beta[2] * c((Tw-1):0, rep(0,25)) + beta[3] * c((Tw-1):0, rep(0,Tp))^2)
    
  }
  
  return(af)
  
}