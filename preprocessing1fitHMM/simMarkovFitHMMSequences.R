lib <- modules::use("R")

library(depmixS4)
library(data.tree)



EstimateModel <- function(dataFrame, outputDir, maxNStates, nPartitions){
  print("estimate model")
  likelihoodsAndTree  <- lib$evalLikelihood$CrossValidationLikelihoodsAndClusteredTree(dataFrame, maxNStates, nPartitions)
  likelihoodsNCluster <- likelihoodsAndTree[[1]]
  plot(1:length(likelihoodsNCluster), likelihoodsNCluster)
  tree <- likelihoodsAndTree[[2]]
  if(likelihoodsNCluster[1] == 0){
    return(list())
  }
  fittedMod <- lib$evalLikelihood$FitMarkovChainFromClusteredTree(dataFrame, tree)
  nStates <- tree$leafCount
  return(list(tree, fittedMod))
}


EstimateDynamicModelParameters <- function(dataFrame, outputDir, maxNStates, nPartitions, i){
  modelList <- EstimateModel(dataFrame, outputDir, maxNStates, nPartitions)
  tree <- modelList[[1]]
  fittedMod <- modelList[[2]]
  nStates <- tree$leafCount
  
  post <- depmixS4::posterior(fittedMod)
  post$executionTime <- dataFrame$outputs
  filename <- paste(outputDir, "executionTimesOccupancyProbabilities", i, ".csv", sep="")
  write.csv(post, file=filename, row.names=FALSE)
  transitionMatrix <- matrix(nrow=nStates, ncol=nStates)
  # get the transition matrix
  for (r in 1:nStates){
    transitionRow = fittedMod@transition[[r]]@parameters$coefficients
    for (c in 1:nStates){
      transitionMatrix[r, c]<- transitionRow[c]
    }
  }
  filename <- paste(outputDir, "transitionMatrix", i, ".txt", sep="")
  write.table(transitionMatrix, file=filename, row.names=FALSE, col.names=FALSE)
  #get the normal distribution mean and stddev
  normalParams <- matrix(nrow=nStates, ncol=2)
  for (r in 1:nStates){
    normalParams[r, 1] <- fittedMod@response[[r]][[1]]@parameters$coefficients[[1]]
    normalParams[r, 2] <- fittedMod@response[[r]][[1]]@parameters$sd
  }
  filename <- paste(outputDir, "normalParams", i, ".txt", sep="")
  write.table(normalParams, file=filename, row.names=FALSE, col.names=FALSE)
  # get the stationary distribution
  # Get the eigenvectors of P, note: R returns right eigenvectors
  r=eigen(transitionMatrix)
  rvec=r$vectors
  # left eigenvectors are the inverse of the right eigenvectors
  lvec=ginv(r$vectors)
  # normalized is the stationary distribution
  # there may be a small imaginary component after inversion, discard that
  pi_eig<-Re(lvec[1,]/sum(lvec[1,]))
  
  filename <- paste(outputDir, "stationaryDistr", i, ".txt", sep="")
  write.table(pi_eig, file=filename, row.names=FALSE, col.names=FALSE)
  
  return()
}

nSequences = 4
preProcessingSequenceLength = 1000
nPartitions <- 4

simSequencesDir = "data/simulatedSequences/"
outputDir <- "data/fittedHMMSequences/"



for(i in 1:4){
  maxNStates <- 4
  filename <- paste(simSequencesDir, "timesStatesSequence", i, ".csv", sep="")
  dataFrame <- read.csv(filename)
  dataFrame <- lib$importData$AdaptDataFrame(dataFrame, 1, TRUE)
  identifyFrame = dataFrame[1:preProcessingSequenceLength,]
  set.seed(3)
  EstimateDynamicModelParameters(identifyFrame, outputDir, maxNStates, nPartitions, i)
}




