adjustment_factor <- function(Tw, Tp, option=1) {
  if (option == 1) {
    
    af = function() rep(1, Tw + Tp)
  
  } else if (option == 2) {
    
    af = function() c(rep(3, Tw-10), rep(2, 10), rep(1, Tp)
    )
    
  } else if (option == 4) {
    
    af = function(beta) exp(beta[1] * ((Tw + Tp - 1):0) + beta[2] * c((Tw-1):0, rep(0,Tp)) + beta[3] * c((Tw-1):0, rep(0,Tp))^2)
    
  } else {
    warning("Invalid option")
    
    return(FALSE)
  }
  
  return(af)
  
}