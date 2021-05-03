library(MASS)
lib <- modules::use("R")


EstimateModel <- function(dataFrame, outputDir, maxNStates, nPartitions, i){
  likelihoodsAndTree  <- lib$evalLikelihood$CrossValidationLikelihoodsAndClusteredTree(dataFrame, maxNStates, nPartitions)
  likelihoodsNCluster <- likelihoodsAndTree[[1]]
  plot(1:length(likelihoodsNCluster), likelihoodsNCluster)
  tree <- likelihoodsAndTree[[2]]
  if(likelihoodsNCluster[1] == 0){
    return(list())
  }
  fittedMod <- lib$evalLikelihood$FitMarkovChainFromClusteredTree(dataFrame, tree)
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
  # get the stationary distribution
  # Get the eigenvectors of P, note: R returns right eigenvectors
  r=eigen(transitionMatrix)
  rvec=r$vectors
  # left eigenvectors are the inverse of the right eigenvectors
  lvec=ginv(r$vectors)
  # normalized is the stationary distribution
  pi_eig<-lvec[1,]/sum(lvec[1,])
  
  filename <- paste(outputDir, "stationaryDistr", i, ".txt", sep="")
  write.table(Re(pi_eig), file=filename, row.names=FALSE, col.names=FALSE)
  
  return(list(tree, fittedMod))
}


posteriorLogLikePart <- function(hyperParams){
  kappas <- hyperParams[,2]
  alphas <- hyperParams[,3]
  betas <- hyperParams[,4]
  res = -lgamma(alphas) + alphas*log(betas) + 0.5*log(kappas)
  return(res)
}

GLRposterior <- function(suffstatsA, suffstatsB, 
                         hyperparamsA, hyperparamsB){
  hyperparams2A <- 
    lib$bayes$updateHyperparamsNormalGamma(suffstatsA, hyperparamsA)
  hyperparams2B <- 
    lib$bayes$updateHyperparamsNormalGamma(suffstatsB, hyperparamsB)
  hyperparamsAB <- 
    lib$bayes$updateHyperparamsNormalGamma(suffstatsA, hyperparamsB)
  suffstatsUnion <- Map("+", suffstatsA, suffstatsB)
  hyperparams2AB <- 
    lib$bayes$updateHyperparamsNormalGamma(suffstatsUnion, hyperparamsAB)
  glr <- -posteriorLogLikePart(hyperparamsA) + posteriorLogLikePart(hyperparams2A) -
    posteriorLogLikePart(hyperparamsB) + posteriorLogLikePart(hyperparams2B) +
    posteriorLogLikePart(hyperparamsAB) - posteriorLogLikePart(hyperparams2AB)
  return(-sum(glr))
}


closestClusterGLRNormalGamma <- function(hyperparamsSegment, hyperparamsClusters,
                                         suffstatsSegment, suffstatsClusters){
  nClusters <- length(hyperparamsClusters)
  if(nClusters == 1){
    return(1)
  }
  minDist <- GLRposterior(suffstatsSegment, suffstatsClusters[[1]],
                          hyperparamsSegment, hyperparamsClusters[[1]])
  minIndex = 1
  for (i in 2:length(hyperparamsClusters)){
    dist <- GLRposterior(suffstatsSegment, suffstatsClusters[[i]], 
                         hyperparamsSegment, hyperparamsClusters[[i]])
    if (dist < minDist) {
      minDist <- dist
      minIndex <- i
    }
  }
  return(minIndex)
}

closestClusterGLRNormalGammaFromInd <- function(hyperparamsClusters, suffstatsClusters,
                                                ind){
  nClusters <- length(hyperparamsClusters)
  if(nClusters == 1){
    return(0, 0)
  }
  startInd = 1
  if (ind == 1){
    startInd = 2
  }
  minDist <- GLRposterior(suffstatsClusters[[ind]], suffstatsClusters[[startInd]],
                          hyperparamsClusters[[ind]], hyperparamsClusters[[startInd]])
  minIndex = startInd
  for (i in startInd:length(hyperparamsClusters)){
    if (i != ind){
      dist <- GLRposterior(suffstatsClusters[[ind]], suffstatsClusters[[i]], 
                           hyperparamsClusters[[ind]], hyperparamsClusters[[i]])
      if (dist < minDist) {
        minDist <- dist
        minIndex <- i
      }
    }
  }
  return(list(minIndex, minDist))
}

closestClustersGLRNormalGammaPreAll <- function(hyperparamsSegment, hyperparamsClusters,
                                                suffstatsSegment, suffstatsClusters,
                                                nPreProcessedClusters){
  nClusters <- length(hyperparamsClusters)
  if(nClusters == 1){
    return(1)
  }
  minDist <- GLRposterior(suffstatsSegment, suffstatsClusters[[1]],
                          hyperparamsSegment, hyperparamsClusters[[1]])
  minIndexPreprocessed = 1
  for (i in 2:nPreProcessedClusters){
    dist <- GLRposterior(suffstatsSegment, suffstatsClusters[[i]], 
                         hyperparamsSegment, hyperparamsClusters[[i]])
    if (dist < minDist) {
      minDist <- dist
      minIndexPreprocessed <- i
    }
  }
  if (nPreProcessedClusters == length(hyperparamsClusters)){
    return(list(minIndexPreprocessed, minIndexPreprocessed))
  }
  minIndexAll <- minIndexPreprocessed
  for (i in (1+nPreProcessedClusters):length(hyperparamsClusters)){
    dist <- GLRposterior(suffstatsSegment, suffstatsClusters[[i]], 
                         hyperparamsSegment, hyperparamsClusters[[i]])
    if (dist < minDist) {
      minDist <- dist
      minIndexAll <- i
    }
  }
  return(list(minIndexPreprocessed, minIndexAll))
}

closestClustersGLRNormalGammaPreAllFromInd <- function(hyperparamsClusters,
                                                       suffstatsClusters, ind,
                                                       nPreProcessedClusters){
  nClusters <- length(hyperparamsClusters)
  if(nClusters == 1){
    return(list(0, 0, 0, 0))
  }
  startInd = 1
  if (ind == 1){
    startInd = 2
  }
  minDistPreprocessed <- GLRposterior(suffstatsClusters[[ind]], suffstatsClusters[[startInd]],
                          hyperparamsClusters[[ind]], hyperparamsClusters[[startInd]])
  minIndexPreprocessed = startInd
  for (i in startInd:nPreProcessedClusters){
    if (i != ind){
      dist <- GLRposterior(suffstatsClusters[[ind]], suffstatsClusters[[i]], 
                           hyperparamsClusters[[ind]], hyperparamsClusters[[i]])
      if (dist < minDistPreprocessed) {
        minDistPreprocessed <- dist
        minIndexPreprocessed <- i
      }
    }
  }
  if (nPreProcessedClusters == length(hyperparamsClusters)){
    return(list(minIndexPreprocessed, minDistPreprocessed, minIndexPreprocessed, minDistPreprocessed))
  }
  minIndexAll <- minIndexPreprocessed
  minDistAll <- minDistPreprocessed
  for (i in (1+nPreProcessedClusters):length(hyperparamsClusters)){
    if (i != ind){
      dist <- GLRposterior(suffstatsClusters[[ind]], suffstatsClusters[[i]], 
                           hyperparamsClusters[[ind]], hyperparamsClusters[[i]])
      if (dist < minDistAll) {
        minDistAll <- dist
        minIndexAll <- i
      }
    }
  }
  return(list(minIndexPreprocessed, minDistPreprocessed, minIndexAll, minDistAll))
}



getLengthOrderedIndices <- function(suffStats){
  lengths <- numeric(length(suffStats))
  for (i in 1:length(suffStats)){
    lengths[i] <- suffStats[[i]][[4]]
  }
  indices <- order(lengths, decreasing=TRUE)
  return(indices)
}


clusterSegmentsBayesian <- function(suffStats, hyperParams, threshold){
  clusterSuffStats <- list()
  lengthOrderedIndices <- getLengthOrderedIndices(suffStats)
  clusterHyperParams <- list()
  clusterSuffStats[[1]] <- suffStats[[lengthOrderedIndices[1]]]
  clusterHyperParams[[1]] <- hyperParams[[lengthOrderedIndices[1]]]
  segmentClusters <- numeric(length(suffStats))
  segmentClusters[lengthOrderedIndices[1]]<-1
  # TODO: cleanup
  for (i in 2:length(suffStats)){
    #    print(clusterSuffStats)
    # find the closest cluster by GLR
    closestIndex <- closestClusterGLRNormalGamma(hyperParams[[lengthOrderedIndices[i]]], clusterHyperParams,
                                                 suffStats[[lengthOrderedIndices[i]]], clusterSuffStats)
    # use GLR to determine whether to merge into this or create a new
    testGLR <- GLRposterior(suffStats[[lengthOrderedIndices[i]]], clusterSuffStats[[closestIndex]], 
                            hyperParams[[lengthOrderedIndices[i]]], clusterHyperParams[[closestIndex]])
    
    if (testGLR < threshold){
      # merge
      segmentClusters[[lengthOrderedIndices[i]]] <- closestIndex
      clusterSuffStats[[closestIndex]] <- Map("+", clusterSuffStats[[closestIndex]], suffStats[[lengthOrderedIndices[i]]])
      clusterHyperParams[[closestIndex]] <- 
        lib$bayes$updateHyperparamsNormalGamma(suffStats[[lengthOrderedIndices[i]]], 
                                                        clusterHyperParams[[closestIndex]])
    }
    else{
      newClusterIndex <- length(clusterSuffStats) + 1
      segmentClusters[[lengthOrderedIndices[i]]] <- newClusterIndex
      clusterSuffStats[[newClusterIndex]] <- suffStats[[lengthOrderedIndices[i]]]
      clusterHyperParams[[newClusterIndex]] <- hyperParams[[lengthOrderedIndices[i]]]
    }
  }
  return(list(segmentClusters, clusterSuffStats, clusterHyperParams))
}


findPotentialSplitSuffStats <-function(suffStatsSteps, clusterHyperParamsLeft, clusterHyperParamsRight, 
                                       suffStatsLeft, suffStatsRight, nSections){
  hyperParamsBegin <- clusterHyperParamsLeft
  beginInd <- 1
  endInd <- nSections - 1
  suffStatsBegin <- suffStatsLeft
  suffStatsEnd <- suffStatsRight
  hyperParamsEnd <- clusterHyperParamsRight
  testBeginSuffStats <- suffStatsSteps[[1]]
  testEndSuffStats <- suffStatsSteps[[nSections]]
  hyperParamsBeginTest <- clusterHyperParamsLeft
  hyperParamsEndTest <- clusterHyperParamsRight
  # TBD: The last section could be split further to find the exact point of change
  while(beginInd < endInd){
    GLRbegin <- GLRposterior(suffStatsBegin, testBeginSuffStats, hyperParamsBegin, hyperParamsBeginTest)
    GLRend <- GLRposterior(suffStatsEnd, testEndSuffStats, hyperParamsEnd, hyperParamsEndTest)
    if (GLRbegin <= GLRend){
      # move interval of possible change forward at beginning
      beginInd <- beginInd + 1
      suffStatsBegin <- Map("+", suffStatsBegin, testBeginSuffStats)
      hyperParamsBegin <- lib$bayes$updateHyperparamsNormalGamma(testBeginSuffStats, hyperParamsBegin)
      testBeginSuffStats <- suffStatsSteps[[beginInd]]
    }
    else{
      # move interval of possible change backwards at end
      endInd <- endInd - 1
      suffStatsEnd <- Map("+", suffStatsEnd, testEndSuffStats)
      hyperParamsEnd <- lib$bayes$updateHyperparamsNormalGamma(testEndSuffStats, hyperParamsEnd)
      testEndSuffStats <- suffStatsSteps[[endInd + 1]]
    }
  }
  return(list(beginInd, hyperParamsBegin, hyperParamsEnd, suffStatsBegin, suffStatsEnd))
}


fbExecTimesOccProbs <- function(df, start, stop, startingProbs, nStates, studentParams, transtionMatrix){
  segment = df[start:stop,]
  forwardResults <- lib$occProbFB$forwardBackwardStudent(segment, nStates, studentParams, 
                                                         startingProbs, transitionMatrix)
  currentExecTimesOccProbs <- data.frame(forwardResults[[1]])
  occProbs <- forwardResults[[2]]
  for (j in 1:nStates){
    colname <- paste("S", j, sep="")
    currentExecTimesOccProbs[colname] <- occProbs[,j]
  }
  currentExecTimesOccProbs["execTime"]<- segment$outputs
  return(currentExecTimesOccProbs)
}


simTestDataDir <- "data/simulatedSequences/"
fittedSequenceDir <- "data/fittedHMMSequences/"
splitIndexDir <- "data/splitIndicesSequencesPreprocess/"
outputDir <- "data/resultsSP/"
nPartitions <- 4
maxNStates <- 4

splittingPointsThresholds = c(20, 20, 20, 20)

for (seq in 1:4){
  mergingThreshold = splittingPointsThresholds[seq]
  print("sequence")
  print(seq)
  filename <- paste(simTestDataDir, "timesStatesSequence", seq, ".csv", sep="")
  dataFrame <- read.csv(filename)
  dataFrame <- lib$importData$AdaptDataFrame(dataFrame, 1, TRUE)
  firstEstimateN <- 1000
  identifyFrame <- dataFrame[1:firstEstimateN,]
  
  set.seed(3)
  # fit with cross validation
  modelList <- EstimateModel(identifyFrame, fittedSequenceDir, maxNStates, nPartitions, seq)
  tree <- modelList[[1]]
  fittedMod <- modelList[[2]]
  etopFilename <- paste(fittedSequenceDir, "executionTimesOccupancyProbabilities", seq, ".csv", sep="")
  execTimesOccProbs <- read.csv(etopFilename)
  weightsFile <- paste(fittedSequenceDir, "stationaryDistr", seq, ".txt", sep="")
  weights <- read.table(weightsFile, header = FALSE, sep = "\n", dec = ".")
  splitIndexFile <- paste(splitIndexDir, "splitIndicesSeq", seq, ".csv", sep="")
  splittingIndices <- read.csv(splitIndexFile)
  nStates <- tree$leafCount
  splittingPoints <- splittingIndices[[1]]
  nSlots <- length(splittingPoints)
  endPreProcessed <- 1000

  normalParams <- matrix(nrow=nStates, ncol=2)
  for (r in 1:nStates){
    normalParams[r, 1] <- fittedMod@response[[r]][[1]]@parameters$coefficients[[1]]
    normalParams[r, 2] <- fittedMod@response[[r]][[1]]@parameters$sd
  }
  
  transitionMatrix <- matrix(nrow=nStates, ncol=nStates)
  # get the transition matrix
  for (r in 1:nStates){
    transitionRow = fittedMod@transition[[r]]@parameters$coefficients
    for (c in 1:nStates){
      transitionMatrix[r, c]<- transitionRow[c]
    }
  }
  
  # cluster segments
  suffStatsList <- list()
  hyperParametersList <- list()
  initHyperParams <- lib$bayes$initializeHyperParametersNormalGamma(normalParams, 20*Re(weights$V1))

  startIndex <- 1
  for(i in 1:nSlots){
    suffStatsList[[i]] <- lib$suffStats$calcSuffStats(execTimesOccProbs, startIndex, 
                                                      splittingPoints[i], nStates)
    startIndex <- splittingPoints[i]+1
  }
  suffStatsList[[nSlots+1]] <- lib$suffStats$calcSuffStats(execTimesOccProbs, startIndex, 
                                                           endPreProcessed, nStates)
  for(i in 1:(nSlots+1)){
    hyperParametersList[[i]] <- 
      lib$bayes$updateHyperparamsNormalGamma(suffStatsList[[i]], 
                                                      initHyperParams)
  }  
  
  clusterSegResult <- clusterSegmentsBayesian(suffStatsList, hyperParametersList, mergingThreshold)
  clusterIndices <- clusterSegResult[[1]]

  outputSplittingPoints <- c(1, splittingPoints)
  
  splittingPointsClustersMat = matrix(c(outputSplittingPoints, clusterIndices), nrow=(length(outputSplittingPoints)))
  filename = paste(outputDir, "clusterIndexSequence", seq, ".csv", sep="")
  write.table(splittingPointsClustersMat, file=filename, row.names=FALSE, col.names=FALSE)
  
  clusterSuffStats <- clusterSegResult[[2]]
  clusterHyperParams <- clusterSegResult[[3]]
  clusterPosteriorStudentParams <- list()

  
  
  for (i in 1:length(clusterSuffStats)){
    clusterPosteriorStudentParams[[i]] <- lib$bayes$studentParams(clusterHyperParams[[i]])
    filename <- paste(outputDir, "studentParams", seq, "_", i, ".txt", sep="")
    write.table(clusterPosteriorStudentParams[[i]], file=filename, row.names=FALSE, col.names=FALSE)
  }

  allStudentParamsList <- c()
  allStudentParamsChanges <- outputSplittingPoints
  for(i in 1:(nSlots)){
    allStudentParamsList <- c(allStudentParamsList, 
                              clusterPosteriorStudentParams[[clusterIndices[[i]]]])
  }  
  
  
  nPreprocessedClusters <- length(clusterSuffStats)

  minLengthSection = 30
  adaptiveProbeChangeThreshold = 0.1 * splittingPointsThresholds[seq]
  adaptiveMergingThreshold = splittingPointsThresholds[seq]
  clusterMergingThresholdPre = 1.5 * splittingPointsThresholds[seq]
  clusterMergingThresholdAll = splittingPointsThresholds[seq]

  slidingWindowStep = 10
  nStepsInWindow = 5
  adaptiveSegmentLength = slidingWindowStep * nStepsInWindow
  
  clusterSearchSegmentLength = 50
  clusterSearchSegmentOffset = 20

  analysisWindowStartIndex = firstEstimateN + 1
  analysisWindowEndIndex = firstEstimateN + adaptiveSegmentLength
  currentClusterIndex = clusterIndices[[length(clusterIndices)]]
  endSequenceIndex = nrow(dataFrame) - clusterSearchSegmentLength + clusterSearchSegmentOffset - slidingWindowStep

  adaptiveStudentParamsList <- c(clusterPosteriorStudentParams[[currentClusterIndex]])
  adaptiveStudentParamsChanges <- c(analysisWindowStartIndex)
  

  initPosteriorStudentParams <- lib$bayes$studentParams(initHyperParams)
  
  startingProbs <- as.numeric(execTimesOccProbs[analysisWindowStartIndex-1, 2:(fittedMod@nstates+1)])
  currentExecTimesOccProbs<- fbExecTimesOccProbs(
    dataFrame, analysisWindowStartIndex, analysisWindowEndIndex, startingProbs, 
    fittedMod@nstates, clusterPosteriorStudentParams[[currentClusterIndex]],
    transtionMatrix)
  suffStatsSteps <- list()
  for(i in 1:nStepsInWindow){
    suffStatsSteps[[i]] <- lib$suffStats$calcSuffStats(currentExecTimesOccProbs, 
                                                       (i-1)*slidingWindowStep + 1, 
                                                       i*slidingWindowStep, fittedMod@nstates)
  }
  
  slidingSuffStats <- suffStatsSteps[[1]]
  for (i in 2:nStepsInWindow){
    slidingSuffStats <- Map("+",slidingSuffStats, suffStatsSteps[[i]])
  }
  
  
  slidingHyperParams <- lib$bayes$updateHyperparamsNormalGamma(slidingSuffStats,
                                                                        initHyperParams)

  while(analysisWindowEndIndex < endSequenceIndex){
    testGLR <- GLRposterior(clusterSuffStats[[currentClusterIndex]], slidingSuffStats, 
                            clusterHyperParams[[currentClusterIndex]], slidingHyperParams)
    changedCluster = FALSE
    # we may have left the current cluster
    if (testGLR > adaptiveProbeChangeThreshold){
      clusterFindSegmentStartIndex <- analysisWindowEndIndex + 1 - clusterSearchSegmentOffset
      clusterFindSegmentEndIndex <- clusterFindSegmentStartIndex + clusterSearchSegmentLength - 1  
      startingProbs <- Re(weights$V1)
      clusterFindETOccProbs<- fbExecTimesOccProbs(
        dataFrame, clusterFindSegmentStartIndex, clusterFindSegmentEndIndex, startingProbs, 
        fittedMod@nstates, initPosteriorStudentParams, transtionMatrix)
      
      clusterFindSuffStats <- lib$suffStats$calcSuffStats(clusterFindETOccProbs, 1, clusterSearchSegmentLength, fittedMod@nstates)
      clusterFindHyperParams <- lib$bayes$updateHyperparamsNormalGamma(clusterFindSuffStats, initHyperParams)
      closestClustersPreAll <- closestClustersGLRNormalGammaPreAll(clusterFindHyperParams, clusterHyperParams,
                                                                   clusterFindSuffStats, clusterSuffStats,
                                                                   nPreprocessedClusters)
      closestClusterPre <- closestClustersPreAll[[1]]

      splitList <- findPotentialSplitSuffStats(suffStatsSteps, clusterHyperParams[[currentClusterIndex]], 
                                               clusterHyperParams[[closestClusterPre]], 
                                               clusterSuffStats[[currentClusterIndex]],
                                               clusterSuffStats[[closestClusterPre]], nStepsInWindow)
      
      splittingInd <- analysisWindowStartIndex + splitList[[1]]*slidingWindowStep
      testEndSuffStats <- suffStatsSteps[[splitList[[1]]]]
      for (i in (splitList[[1]] + 1):nStepsInWindow){
        testEndSuffStats <- Map("+", testEndSuffStats, suffStatsSteps[[i]])
      }
      
      testNewHyperParamsPre <- lib$bayes$updateHyperparamsNormalGamma(testEndSuffStats, clusterHyperParams[[closestClusterPre]])
      if (closestClusterPre != currentClusterIndex){ 
        #change current cluster
        currentClusterIndex <- closestClusterPre
        # add splitting point
        splittingPoints[[length(splittingPoints)+1]] <- splittingInd
        clusterIndices[[length(splittingPoints)+1]] <- currentClusterIndex
        adaptiveStudentParamsList <- c(adaptiveStudentParamsList, clusterPosteriorStudentParams[[currentClusterIndex]])
        adaptiveStudentParamsChanges <- c(adaptiveStudentParamsChanges, splittingInd)
        changedCluster = TRUE
      } 
    }
    if (changedCluster){
      # move until end of analysis window when we have had a split
      startingProbs <- Re(weights$V1)
      analysisWindowStartIndex <- analysisWindowEndIndex+1
      analysisWindowEndIndex <- analysisWindowEndIndex+slidingWindowStep*nStepsInWindow
      if (analysisWindowEndIndex < endSequenceIndex){
        currentExecTimesOccProbs<- fbExecTimesOccProbs(
          dataFrame, analysisWindowStartIndex, analysisWindowEndIndex, startingProbs, 
          fittedMod@nstates, clusterPosteriorStudentParams[[currentClusterIndex]],
          transtionMatrix)
        suffStatsSteps <- list()
        for(i in 1:nStepsInWindow){
          suffStatsSteps[[i]] <- lib$suffStats$calcSuffStats(currentExecTimesOccProbs, 
                                                             (i-1)*slidingWindowStep + 1, 
                                                             i*slidingWindowStep, fittedMod@nstates)
        }
        slidingSuffStats <- suffStatsSteps[[1]]
        for (i in 2:nStepsInWindow){
          slidingSuffStats <- Map("+",slidingSuffStats, suffStatsSteps[[i]])
        }
        slidingHyperParams <- lib$bayes$updateHyperparamsNormalGamma(slidingSuffStats,
                                                                              initHyperParams)
      }
    }
    else{
      # step
      startingProbs <- as.numeric(currentExecTimesOccProbs[(slidingWindowStep*nStepsInWindow), 2:(fittedMod@nstates+1)])
      stepOccProbs<- fbExecTimesOccProbs(
        dataFrame, analysisWindowEndIndex+1, analysisWindowEndIndex+slidingWindowStep, startingProbs, 
        fittedMod@nstates, clusterPosteriorStudentParams[[currentClusterIndex]],
        transtionMatrix)
      
      stepSuffStats <- lib$suffStats$calcSuffStats(stepOccProbs, 1, slidingWindowStep, fittedMod@nstates)
      
      currentExecTimesOccProbs[1:(slidingWindowStep*(nStepsInWindow - 1)),] <- 
        currentExecTimesOccProbs[(slidingWindowStep + 1):(slidingWindowStep*(nStepsInWindow)),]
      currentExecTimesOccProbs[(slidingWindowStep*(nStepsInWindow-1)+1):(slidingWindowStep*nStepsInWindow),] <- 
        stepOccProbs
      stepSuffStats <- lib$suffStats$calcSuffStats(stepOccProbs, 1, slidingWindowStep, fittedMod@nstates)
      
      
      removeSuffStatsStep <- suffStatsSteps[[1]]
      suffStatsSteps[1:(nStepsInWindow - 1)] <- suffStatsSteps[2:nStepsInWindow]
      suffStatsSteps[[nStepsInWindow]] <- stepSuffStats
      
      suffStatsDiff <- Map("-",stepSuffStats, removeSuffStatsStep)
      slidingSuffStats <- Map("+",slidingSuffStats, suffStatsDiff)
      slidingHyperParams <- lib$bayes$updateHyperparamsNormalGamma(stepSuffStats,
                                                                            slidingHyperParams)
      slidingHyperParams <- lib$bayes$updateHyperparamsNormalGammaNeg(removeSuffStatsStep,
                                                                               slidingHyperParams)
      analysisWindowStartIndex <- analysisWindowStartIndex + slidingWindowStep
      analysisWindowEndIndex <- analysisWindowEndIndex + slidingWindowStep
    }
  }
  filename = paste(outputDir, "adaptiveStudentParams", seq, ".csv", sep="")
  write.table(as.data.frame(adaptiveStudentParamsList, stringsAsFactors = FALSE),file=filename, quote=F,sep=",",row.names=F, col.names=F)
  
  filename = paste(outputDir, "adaptiveStudentParamsChangeIndices", seq, ".csv", sep="")
  write.table(as.data.frame(adaptiveStudentParamsChanges, stringsAsFactors = FALSE),file=filename, quote=F,sep=",",row.names=F, col.names=F)
  
  outputSplittingPoints <- c(1, splittingPoints)
  allStudentParamsList <- c(allStudentParamsList, adaptiveStudentParamsList)
  allStudentParamsChanges <- c(allStudentParamsChanges, adaptiveStudentParamsChanges[-1])

  splittingPointsClustersMat = matrix(c(outputSplittingPoints, clusterIndices), nrow=(length(outputSplittingPoints)))
  filename = paste(outputDir, "clusterIndexSequenceFull", seq, ".csv", sep="")
  write.table(splittingPointsClustersMat, file=filename, row.names=FALSE, col.names=FALSE)
  
  for (i in 1:length(clusterSuffStats)){
    clusterPosteriorStudentParams[[i]] <- lib$bayes$studentParams(clusterHyperParams[[i]])
    filename <- paste(outputDir, "postAdaptive/studentParams", seq, "_", i, ".txt", sep="")
    write.table(clusterPosteriorStudentParams[[i]], file=filename, row.names=FALSE, col.names=FALSE)
  }
}






