

calc_expected_CV <- function (strata) {
  
  n_spp <- length(grep(x = names(strata), pattern = "M"))
  n_h <- strata$Allocation
  N_h <- strata$Population
  
  cv <- vector(length = n_spp)
  names(cv) <- paste0("Y", 1:n_spp)
  
  for (ispp in 1:n_spp) {
    S_h <- strata[, paste0("S", ispp)]
    M_h <- strata[, paste0("M", ispp)]
    
    Y_h <- N_h * M_h
    Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
    CV <- sqrt(sum(Var_h))/sum(Y_h)
    
    cv[ispp] <- CV
  }
  
  cv <- round(cv, 3)
  return(cv)
}