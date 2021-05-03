library(ggplot2)
library(MASS)



nSequences = 4
sequenceLength = 3000
# a new cluster can be added after this period
preProcessingSequenceLength = 1000
nStates = 3 # don't change, transition matrix generation has this fixed
nClustersPreprocesing = 4
nAdditionalClusters = 1
segmentLengthMin = 50
segmentLengthMax = 300

statesMeanLoRange = c(25, 60, 95)   
statesMeanHiRange = c(50, 85, 120)
statesStdDevLoRange = c(2, 2, 2)
statesStdDevHiRange = c(6, 6, 6)

outputDir = "data/simulatedSequences/"

set.seed(3)

orderedMeans =  list()
for (i in 1:nStates){
  orderedMeansVec = sort(runif(nClustersPreprocesing + nAdditionalClusters, 
                               min = statesMeanLoRange[i],
                               max = statesMeanHiRange[i]))
  orderedMeans[[i]] = orderedMeansVec
}

clusterMapping <- sample.int(5)
print(clusterMapping)
print(orderedMeans)


for (i in 1:nSequences){
  means = matrix(nrow = nClustersPreprocesing + nAdditionalClusters,
                 ncol = nStates)
  stdDevs = matrix(nrow = nClustersPreprocesing + nAdditionalClusters,
                   ncol = nStates)
  transitionMatrix = matrix(nrow=3, ncol=3)
  transitionMatrix[,1] = runif(3, min=0.1, max=0.8)

  for (state in 1:nStates){
    for (c in 1:(nClustersPreprocesing + nAdditionalClusters)){
      means[c,state] = orderedMeans[[state]][clusterMapping[c]]
    }
    stdDevs[,state] = runif(nClustersPreprocesing + nAdditionalClusters, 
                            min = statesStdDevLoRange[state],
                            max = statesStdDevHiRange[state])
    transitionMatrix[state,2] = runif(1, min=0.1, max=(1-transitionMatrix[state,1]-0.1))
    transitionMatrix[state,3] = 1 - sum(transitionMatrix[state,1:2])
  }
  filename = paste(outputDir, "meansSequence", i, ".csv", sep="")
  write.table(means, file=filename, row.names=FALSE, col.names=FALSE)
  filename = paste(outputDir, "stdDevsSequence", i, ".csv", sep="")
  write.table(stdDevs, file=filename, row.names=FALSE, col.names=FALSE)
  filename = paste(outputDir, "transitionMatrix", i, ".csv", sep="")
  write.table(transitionMatrix, file=filename, row.names=FALSE, col.names=FALSE)
  r=eigen(transitionMatrix)
  rvec=r$vectors
  # left eigenvectors are the inverse of the right eigenvectors
  lvec=ginv(r$vectors)
  # normalized is the stationary distribution
  pi_eig<-lvec[1,]/sum(lvec[1,])
  
  filename <- paste(outputDir,"stationaryDistr", i, ".csv", sep="")
  write.table(Re(pi_eig), file=filename, row.names=FALSE, col.names=FALSE)
  
  splittingPoints = c(1)
  preProcessingClusters = 1:nClustersPreprocesing
  currentCluster = sample(preProcessingClusters, 1)
  clusters = c(currentCluster)
  currentPoint = 1
  while(currentPoint < preProcessingSequenceLength){
    currentPoint = currentPoint + sample(segmentLengthMin:segmentLengthMax, 1)
    splittingPoints = c(splittingPoints, currentPoint)
    switchToClustersIndices = which(preProcessingClusters != currentCluster)
    currentCluster = sample(preProcessingClusters[switchToClustersIndices], 1)
    clusters = c(clusters, currentCluster)
  }
  allClusters = 1:(nClustersPreprocesing + nAdditionalClusters)
  while(currentPoint < sequenceLength - segmentLengthMax){
    currentPoint = currentPoint + sample(segmentLengthMin:segmentLengthMax, 1)
    splittingPoints = c(splittingPoints, currentPoint)
    switchToClustersIndices = which(allClusters != currentCluster)
    currentCluster = sample(allClusters[switchToClustersIndices], 1)
    clusters = c(clusters, currentCluster)
  }
  splittingPointsClustersMat = matrix(c(splittingPoints, clusters), nrow=(length(splittingPoints)))
  filename = paste(outputDir, "splittingPointsClustersSequence", i, ".csv", sep="")
  write.table(splittingPointsClustersMat, file=filename, row.names=FALSE, col.names=FALSE)

  splittingPoints = c(splittingPoints, sequenceLength)
  executionTime = numeric(sequenceLength)
  stateSeq = numeric(sequenceLength)
  currentState = 1

  for(segment in 1:length(clusters)){
    for(sample in splittingPoints[segment]:splittingPoints[segment+1]){
      stateProbs = transitionMatrix[currentState,]
      currentState = sample(1:3, size=1, prob = stateProbs)
      stateSeq[sample] = currentState
      executionTime[sample] = rnorm(1, mean= means[clusters[segment], currentState], 
                              sd = stdDevs[clusters[segment], currentState])
    }     
  }
  
  df = data.frame(executionTime)
  chainLength <- length(executionTime)
  df$index = seq(1, chainLength)
  df$states = stateSeq
  filename = paste(outputDir, "timesStatesSequence", i, ".csv", sep="")
  write.csv(df, filename, row.names = FALSE)

  gg <- ggplot(df, aes(y = executionTime, x = index)) + geom_point() 
    
  filename = paste(outputDir, "timesStatesSequenceImage", i, ".png", sep="")  
  ggsave(filename)
  
}

