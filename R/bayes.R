
studentParams <- function(hyperParams){
  studentParams <- matrix(nrow=nrow(hyperParams), ncol=3)
  means <- hyperParams[,1]
  scaling <- sqrt(hyperParams[,4]*(hyperParams[,2] + 1)/(hyperParams[,2]*hyperParams[,3]))
  dof <- hyperParams[,2]
  studentParams[,1] <- means
  studentParams[,2] <- scaling
  studentParams[,3] <- dof
  return(studentParams)
}



initializeHyperParametersNormalGamma <- function(normalParams, nPseudoObs){
  nStates <- nrow(normalParams)
  hyperParams <- matrix(nrow=nStates, ncol=4)
  # mu - mean
  hyperParams[,1] <- normalParams[,1]
  # nu - number of pseudo obs for mean estimate
  hyperParams[,2] <- nPseudoObs
  # alpha - number of pseudo obs for sum of squared dev/2
  hyperParams[,3] <- nPseudoObs/2
  # beta - sum of squared deviations/2
  hyperParams[,4] <- hyperParams[,3] * normalParams[,2]^2
  return(hyperParams)
}


updateHyperparamsNormalGamma <- function(suffStats, hyperParams, threshold = 0.01){
  a0 <- suffStats[[1]]
  a1 <- suffStats[[2]]
  a2 <- suffStats[[3]]
  valid <- a0 > threshold
  hyperParams_out <- matrix(nrow=nrow(hyperParams), ncol=4)
  hyperParams_out[,] <- hyperParams[,]
  mu_in <- hyperParams[,1]
  nu_in <- hyperParams[,2]
  # alpha is redundant - nu/2
  alpha_in <- hyperParams[,3]
  beta_in <- hyperParams[,4]
  # mu update
  hyperParams_out[valid,1] <- (nu_in[valid]*mu_in[valid] + a1[valid])/(nu_in[valid] + a0[valid])
  # nu update
  hyperParams_out[valid, 2] <- nu_in[valid] + a0[valid]
  # alpha is redundant - nu/2
  hyperParams_out[valid,3] <- alpha_in[valid] + a0[valid]/2
  # beta update
  hyperParams_out[valid, 4] <- beta_in[valid] + 
    0.5*(a2[valid] - a1[valid]^2/a0[valid] + 
           nu_in[valid] * a0[valid] * (a1[valid]/a0[valid] - mu_in[valid])^2/ 
           (nu_in[valid] + a0[valid]))
  return(hyperParams_out)
}

updateHyperparamsNormalGammaNeg <- function(suffStats, hyperParams, threshold = 0.01){
  a0 <- suffStats[[1]]
  a1 <- suffStats[[2]]
  a2 <- suffStats[[3]]
  valid <- a0 > threshold
  hyperParams_out <- matrix(nrow=nrow(hyperParams), ncol=4)
  hyperParams_out[,] <- hyperParams[,]
  mu_in <- hyperParams[,1]
  nu_in <- hyperParams[,2]
  # alpha is redundant - nu/2
  alpha_in <- hyperParams[,3]
  beta_in <- hyperParams[,4]
  # mu update
  hyperParams_out[valid,1] <- (nu_in[valid]*mu_in[valid] - a1[valid])/(nu_in[valid] - a0[valid])
  # nu update
  hyperParams_out[valid, 2] <- nu_in[valid] - a0[valid]
  # alpha is redundant - nu/2
  hyperParams_out[valid,3] <- alpha_in[valid] - a0[valid]/2
  # beta update
  hyperParams_out[valid, 4] <- beta_in[valid] - 
    0.5*(a2[valid] - a1[valid]^2/a0[valid] + 
           nu_in[valid] * a0[valid] * (a1[valid]/a0[valid] - mu_in[valid])^2/ 
           (nu_in[valid] - a0[valid]))
  return(hyperParams_out)
}


