
export("forwardBackwardStudent")
forwardBackwardStudent<- function(dataFrame, nStates, studentParams, initProbs, transitionMatrix){
  data <- dataFrame$outputs
  alphas <- matrix(nrow=length(data), ncol=nStates)
  scaling <- matrix(nrow=length(data), ncol=1)
  betas <- matrix(nrow=length(data), ncol=nStates)
  gammas <- matrix(nrow=length(data), ncol=nStates)
  studentStateProbs <-  matrix(nrow=length(data), ncol=nStates)
  for (i in 1:length(data)){
    for (j in 1:nStates){
      studentStateProbs[i,j] <- stats::dt( 
        (data[i] - studentParams[j,1])/ studentParams[j,2], df=studentParams[j,3])
    }
  }
  for (i in 1:nStates){
    alphas[1,i] <- initProbs[i]* studentStateProbs[1,i]
    betas[length(data),i] <- 1
  }
  scaling[1,1] <- 1/sum(alphas[1,])
  alphas[1,] <- alphas[1,]*scaling[1,1]
  
  for (i in 2:length(data)){
    for (j in 1:nStates){
      probj = sum(transitionMatrix[,j]*alphas[i-1,])
      alphas[i,j] <- probj* studentStateProbs[i,j]
    }
    scaling[i,1] <- 1/sum(alphas[i,])
    alphas[i,] <- alphas[i,]*scaling[i,1]
  }
  gammas[length(data),] = alphas[length(data),]
  for (i in (length(data)-1):1){
    for (j in 1:nStates){
      betas[i,j] = sum(transitionMatrix[j,]*betas[i+1,]*studentStateProbs[i+1,])
    }
    betas[i,] = betas[i,]*scaling[i,1]
    gammas[i,] = betas[i,]*alphas[i,]
    gammas[i, ] = gammas[i,]/sum(gammas[i,])
  }
  return(list(scaling, gammas))
}




