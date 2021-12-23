pension_fund <- function(x, E, nSim, T=600, afx=1, gamma=5, p=0.2) {
  #' pension_fund() simulates a pension fund with portfolio weights x, 
  #' the function returns the certainty equivalent of consumption based on
  #' T years and nSim simulations.
  #'
  #' @param x: weight of portfolio
  #' @param E: list with economic data
  #' @param nSim: number of simulations
  #' @param afx: adjustment factor option
  #'
  #' @return: certainty equivalent of consumption (CEC)
  #'
  #' @examples
  
  Tw = 40
  Tp = 25
  
  # Q-matrix
  Q = matrix(0, Tw + Tp, Tw + Tp)
  for (i in 1:(Tw+Tp)) {
    for (j in 1:(Tw+Tp)) {
      if (((i + j) > (Tw + 1)) & ((i + j) < (Tw + Tp + 2))) {
        Q[i,j] = 1
      }
    }
  }
  
  # Create B-matrix, everyone receives 40% of their salary as pension
  B = array(0, c(Tw + Tp, Tw + Tp, T+1))
  disc = .4 * c(seq(1/Tw, 1, by=1/Tw), rep(1, Tp))
  B[,,1] = sweep(Q, 2, E$w[1] * disc, "*")
  
  A           = matrix(0, T+1, 1)
  L           = matrix(0, T+1, 1)
  discB       = array(0, c(Tw+Tp, Tw+Tp, T+1))
  Ink         = matrix(0, T, 1)
  Uit         = matrix(0, T, 1)
  Obl         = matrix(0, T, 1)
  CR          = matrix(1, T, 1)
  v           = array(0, c(Tw+Tp, Tw+Tp, T+1))
  ind         = matrix(0, T, 1)
  ksi         = matrix(0, T, 1)
  U           = matrix(0, T-Tp+1, nSim)
  CRstrike    = 0
  afag        = matrix(0, 1, Tw+Tp)
  
  alpha = matrix(1, Tw+Tp, 1)
  for (i in 1:9) {
    alpha[i] = 1/10
  }
  
  af = adjustment_factor(Tw, Tp, afx)
  af = if(afx==4) af(beta) else af()
  
  # Needed for utility
  u       = function(x, gamma=5) x^(1-gamma) / (1-gamma)
  
  # Rho option 1
  # rho     = 1 / (1 + mean(E$r))
  
  # Rho option 2
  Xagg    = t(sapply(1:2, function(i) colSums(E$X[i,,]))) / nSim
  rho     = 1 / (1 + E$delta_r[1] + mean(E$delta_r[2:3] %*% Xagg))
    
  B[,,2] = B[,,1]
  
  for (s in 1:nSim) {
    print(paste0("Running simulation: ", s,"/",nSim))
    discB[,,1] = E$P[,s,1] %*% B[,,1]
    discB[,,2] = E$P[,s,2] %*% B[,,2]
    
    # TODO: build adjustment factor for adjust.. achievable pension..  
    if (afx == 3) {
      
    }
    
    L[1] = sum(discB[,,1])
    A[1] = L[1]
    
    for (year in 2:T) {
      # Update assets and liabilities, and compute new coverage ratio
      B[,,year]     = B[,,(year-1)]
      discB[,,year] = B[,,year] * E$P[,s,year]
      L[year]       = sum(discB[,,year])
      A[year]       = A[year-1] * ( x[1] * (1 + E$dS_[year,s]) + (1 - x[1] - x[2]) * (1 + E$r[year,s]) ) + x[2] * CR[year-1] * L[year]
      CR[year]      = A[year] / L[year]
      
      # Add penalty if coverage ratio is under 100%
      if (CR[year] < 1) {
        CRstrike = CRstrike + 1
      } else {
        CRstrike = 0
      }
      
      if (CRstrike == 5) {
        ind[year] = CR[year] - 1
        # ksi[year]  = ind[year] * L[year] / sum(sweep(sweep(B[,,year], 2, alpha * af, "*"), 1, E$P[,s,year], "*"))
        ksi[year]  = ind[year] * L[year] / sum(E$P[,s,year] * t(t(B[,,year]) * c(alpha*af)))
        v[,,year]  = ksi[year] * alpha %*% t(af)
      } else if ( CR[year] > 0.9 & CR[year] < 1.2) {
        ind[year] = (A[year] - L[year]) / (9 * A[year] + L[year])
        ksi[year] = ind[year] * L[year] / sum(E$P[,s,year] * t(t(B[,,year]) * af))
        v[,,year] = ksi[year] * rep(1, Tw + Tp) %*% t(af)
      } else if ( CR[year] > 1.2) {
        ind[year] = (A[year] - L[year]) / (4 * A[year] + L[year])
        ksi[year] = ind[year] * L[year] / sum(E$P[,s,year] * t(t(B[,,year]) * af))
        v[,,year] = ksi[year] * rep(1, Tw + Tp) %*% t(af)
      } else {
        ind[year] = CR[year] / 0.9 - 1
        ksi[year] = ind[year] * L[year] / sum(E$P[,s,year] * t(t(B[,,year]) * c(alpha*af)))
        v[,,year] = ksi[year] * alpha %*% af
      }
      
      Ink[year] = Tw * p * E$w[year,s];
      B[,,year] = (1 + v[,,year]) * B[,,year]
      
      # Nonnegativity restriction for pensions
      B[,,year][B[,,year] < 0] = 0
      # B[,,year] = (abs(B[,,year]) + B[,,year])/2
      
      discB[,,year] = E$P[,s,year] * B[,,year]
      L[year]       = sum(discB[,,year])
      CR[year]      = A[year] / L[year]
      Uit[year]     = sum(B[1, (Tw+1):(Tw+Tp), year])
      
      A[year] = A[year] - Uit[year] + Ink[year]
      B[1:(Tw+Tp-1), 2:(Tw+Tp), year] = B[2:(Tw+Tp), 1:(Tw+Tp-1), year]
      B[(Tw+Tp),,year] = 0
      B[,1,year] = 0
      
      # Update outstanding payoffs with incoming wage
      B[,,year] = B[,,year] + t(t(Q) * ((p*E$w[year,s]) / colSums(E$P2[,s,year] * Q)))
      
      discB[,,year] = B[,,year] * E$P2[,s,year]
      L[year] = sum(discB[,,year])
      CR[year] = A[year] / L[year]   
      
      # TODO
      if (afx == 3) {
        
      }
      
    } # end time loop
    
    # Compute utility
    rho_ = 1 / (mean(E$r[,s]) + 1)
    for (j in (-Tw:(T-Tw-Tp))) {
      U[j+Tw+1,s] = 0
      for (t in (j+Tw+1):(j+Tw+Tp)) {
        U[j+Tw+1, s] = U[j+Tw+1, s] + rho_^(t-j-Tw-1) * u(B[1,t-j,t] / E$Pi[t,s])
      }
    }
  } # end simulation loop
  
  welfare = rowSums(U)
  SW = sum(rho^(100:length(welfare)) * welfare[100:length(welfare)])
  CEC = -((SW*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1))
  
}