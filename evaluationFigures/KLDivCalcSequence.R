library(ggplot2)
library(LaplacesDemon)


calcKLDiv <- function(studentMeans, studentScales, studentDofs, studentProbs, 
                          normalMeans, normalStddevs, normProbs, n){
  estProbs = numeric(length = n)
  x = runif(n, min=0, max=150)
  for (i in 1:nrow(studentProbs)){
    mean <- studentMeans[i]
    scale <- studentScales[i]
    dof <- studentDofs[i]
    estProbs = estProbs + studentProbs[i,1]*
      stats::dt((x - mean)/scale, df=dof)
  }
  trueProbs = numeric(length = n)
  for (i in 1:nrow(normProbs)){
    trueProbs = trueProbs + normProbs[i,1]*
      stats::dnorm(x, mean = normalMeans[i], sd= normalStddevs[i])
  }
  res = KLD(trueProbs, estProbs)
  return(res$sum.KLD.px.py)
}


getAllPreprocessingChangePoints <- function(trueInputDir, estInputDir, indexIn){
  trueSplittingPointsFile = paste(trueInputDir, "splittingPointsClustersSequence", indexIn, ".csv", sep="")
  trueSplittingPoints = read.table(trueSplittingPointsFile)
  trueSplitIndices = trueSplittingPoints[,1]
  preIndices = which(trueSplitIndices <= 1000)
  changingPointsClusters = paste(estInputDir, "clusterIndexSequence", indexIn, ".csv", sep="")
  pointsClusters = read.table(changingPointsClusters)
  allChangePoints = c(trueSplitIndices[preIndices], pointsClusters[,1])
  s = sort(allChangePoints)
  return(s)
}

getAllAdaptiveChangePoints <- function(trueInputDir, estInputDir, indexIn){
  trueSplittingPointsFile = paste(trueInputDir, "splittingPointsClustersSequence", indexIn, ".csv", sep="")
  trueSplittingPoints = read.table(trueSplittingPointsFile)
  trueSplitIndices = trueSplittingPoints[,1]
  adaptiveStudentChangesFile = paste(estInputDir, "adaptiveStudentParamsChangeIndices", indexIn, ".csv", sep="")
  adaptiveChangePoints = read.table(adaptiveStudentChangesFile, header=FALSE, stringsAsFactors= FALSE, sep=",")
  trueAdaptiveIndices = which(trueSplitIndices > 1000)
  allChangePoints = c(trueSplitIndices[trueAdaptiveIndices], adaptiveChangePoints[,1])
  s = sort(allChangePoints)
  return(s)
}


getNormalParams <- function(trueInputDir, indexIn){
  trueSplittingPointsFile = paste(trueInputDir, "splittingPointsClustersSequence", indexIn, ".csv", sep="")
  trueSplittingPoints = read.table(trueSplittingPointsFile)
  trueMeansFile = paste(trueInputDir, "meansSequence", indexIn, ".csv", sep="")
  trueMeans = read.table(trueMeansFile)
  trueStddevsFile = paste(trueInputDir, "stdDevsSequence", indexIn, ".csv", sep="")
  trueStdDevs = read.table(trueStddevsFile)
  trueSplittingPointIndex = 2
  meansMatrix = matrix(nrow = 3000, ncol=3)
  stddevsMatrix = matrix(nrow = 3000, ncol=3)
  clustersMatrix = matrix(nrow = 3000, ncol=1)
  startIndex = 1
  for (ind in trueSplittingPointIndex:nrow(trueSplittingPoints)){
    stopIndex = trueSplittingPoints[ind,1]-1
    for (i in startIndex:stopIndex){
      for (j in 1:3){
        meansMatrix[i,j] = trueMeans[trueSplittingPoints[ind-1,2],j]
        stddevsMatrix[i,j] = trueStdDevs[trueSplittingPoints[ind-1,2],j]
      }
      clustersMatrix[i,1] = trueSplittingPoints[ind-1,2]
    }
    startIndex = stopIndex+1
  }
  for (i in startIndex:3000){
    for (j in 1:3){
      meansMatrix[i,j] = trueMeans[trueSplittingPoints[length(trueSplittingPoints),2],j]
      stddevsMatrix[i,j] = trueStdDevs[trueSplittingPoints[length(trueSplittingPoints),2],j]
    }
    clustersMatrix[i,1] = trueSplittingPoints[length(trueSplittingPoints),2]
    
  }
  return(list(meansMatrix, stddevsMatrix, clustersMatrix))
}


addAdaptiveStudentParams <- function(adaptiveInputDir, indexIn, preprocessingResult){
  adaptiveStudentChangesFile = paste(adaptiveInputDir, "adaptiveStudentParamsChangeIndices", indexIn, ".csv", sep="")
  adaptiveChangePoints = read.table(adaptiveStudentChangesFile, header=FALSE, stringsAsFactors= FALSE, sep=",")
  adaptiveStudentParamsFile = paste(adaptiveInputDir, "adaptiveStudentParams", indexIn, ".csv", sep="")
  adaptiveStudentParams = read.table(adaptiveStudentParamsFile, header=FALSE, stringsAsFactors= FALSE, sep=",")
  nStates = nrow(adaptiveStudentParams)/(3*nrow(adaptiveChangePoints))
  studentMeans = preprocessingResult[[1]]
  studentScales = preprocessingResult[[2]]
  studentDofs = preprocessingResult[[3]]
  firstChangePoint = adaptiveChangePoints[1, 1]
  for (i in 1001:(firstChangePoint-1)){
    studentMeans[i,] = studentMeans[1000,]
    studentScales[i,] = studentScales[1000,]
    studentDofs[i,] = studentDofs[1000,]
  }
  startIndex = firstChangePoint
  for(row in 1:(nrow(adaptiveChangePoints)-1)){
    startReadIndexParams = (row-1)*(3*nStates)
    meansIndex = startReadIndexParams
    stopIndex = adaptiveChangePoints[row+1, 1]-1
    for (state in 1:nStates){
      meansIndex = meansIndex + 1
      scaleIndex = meansIndex + nStates
      dofIndex = scaleIndex + nStates
      mean = adaptiveStudentParams[meansIndex,1]
      sc = adaptiveStudentParams[scaleIndex,1]
      dof = adaptiveStudentParams[dofIndex,1]
      studentMeans[startIndex:stopIndex, state] = mean
      studentScales[startIndex:stopIndex, state] = sc
      studentDofs[startIndex:stopIndex, state] = dof
    }
    startIndex = stopIndex + 1
  }
  startReadIndexParams = (nrow(adaptiveChangePoints)-1)*(3*nStates)
  meansIndex = startReadIndexParams
  stopIndex = 3000
  for (state in 1:nStates){
    meansIndex = meansIndex + 1
    scaleIndex = meansIndex + nStates
    dofIndex = scaleIndex + nStates
    mean = adaptiveStudentParams[meansIndex,1]
    sc = adaptiveStudentParams[scaleIndex,1]
    dof = adaptiveStudentParams[dofIndex,1]
    studentMeans[startIndex:stopIndex, state] = mean
    studentScales[startIndex:stopIndex, state] = sc
    studentDofs[startIndex:stopIndex, state] = dof
  }
  return(list(studentMeans, studentScales, studentDofs))
}

getPreprocessingStudentParams <- function(preprocessingInputDir, indexIn){
  changingPointsClusters = paste(preprocessingInputDir, "clusterIndexSequence", indexIn, ".csv", sep="")
  pointsClusters = read.table(changingPointsClusters)
  nClusters = max(pointsClusters[,2])
  studentParams = list()
  for (cluster in 1:nClusters){
    studentParamsFileName = paste(preprocessingInputDir, "studentParams", indexIn, "_", cluster, ".txt", sep="")
    studentParams[[cluster]] = read.table(studentParamsFileName)
  }
  nStates = nrow(studentParams[[1]])
  studentMeans = matrix(nrow=3000, ncol=nStates)
  studentScales = matrix(nrow=3000, ncol=nStates)
  studentDofs = matrix(nrow=3000, ncol=nStates)
  startIndex = pointsClusters[1,1]
  for(row in 1:(nrow(pointsClusters)-1)){
    stopIndex = pointsClusters[row+1, 1]-1
    studentParamsCluster = studentParams[[pointsClusters[row,2]]]
    for (state in 1:nStates){
      mean = studentParamsCluster[state, 1]
      scale = studentParamsCluster[state, 2]
      dof = studentParamsCluster[state, 3]
      studentMeans[startIndex:stopIndex, state] = mean
      studentScales[startIndex:stopIndex, state] = scale
      studentDofs[startIndex:stopIndex, state] = dof
    }
    startIndex = stopIndex + 1
  }
  stopIndex = 1000
  studentParamsCluster = studentParams[[pointsClusters[nrow(pointsClusters),2]]]
  for (state in 1:nStates){
    mean = studentParamsCluster[state, 1]
    scale = studentParamsCluster[state, 2]
    dof = studentParamsCluster[state, 3]
    studentMeans[startIndex:stopIndex, state] = mean
    studentScales[startIndex:stopIndex, state] = scale
    studentDofs[startIndex:stopIndex, state] = dof
  }
  return(list(studentMeans, studentScales, studentDofs))
}

getAllStudentParams<-function(inputDir, indexIn){
  result = getPreprocessingStudentParams(inputDir, indexIn)
  return(addAdaptiveStudentParams(inputDir, indexIn, result))
}


outputFormattedKLDivSequence <- function(KLDivMatrixList, seqIndex, outputDir){
  filename = paste(outputDir, "KLDiv_", seqIndex, ".txt", sep="")
  cat("", file = filename, append = FALSE)
  for (cluster in 1:6){
    for (process in 1:4){
      klDivMean = KLDivMatrixList[[process]][cluster, seqIndex]
      cat(sprintf(' & %0.3f', klDivMean), file = filename, append = TRUE)
    }
    cat(sprintf('\n'), file = filename, append = TRUE)
  }
}


linesMatrixAdaptive <- function(trueInputDir, fittedModelDir, adaptiveInputDir, 
                                outputDir, outputFileString){
  clustersKLDivMat = matrix(0, nrow = 6, ncol=4)
  for (i in 1:4){
    trueStateProbsFile = paste(trueInputDir, "stationaryDistr", i, ".csv", sep="")
    trueStateProbs = read.table(trueStateProbsFile)
    estStateProbsFile = paste(fittedModelDir, "stationaryDistr", i, ".txt", sep="")
    estStateProbs = read.table(estStateProbsFile)
    normalParams = getNormalParams(trueInputDir, i)
    studentParams = getAllStudentParams(adaptiveInputDir, i)
    changePoints = getAllAdaptiveChangePoints(trueInputDir, adaptiveInputDir, i)
    preLines = matrix(nrow = 2*length(changePoints), ncol= 2)
    nPerCluster = numeric(5)
    startingPoint = 1001
    for (j in 1:length(changePoints)){
      stoppingPoint = changePoints[[j]]-1
      preLines[(j-1)*2-1,1] = startingPoint
      kl = calcKLDiv(studentParams[[1]][startingPoint,], studentParams[[2]][startingPoint,], 
                     studentParams[[3]][startingPoint,], estStateProbs, 
                     normalParams[[1]][startingPoint,], normalParams[[2]][startingPoint,], 
                     trueStateProbs, 10000)
      preLines[(j-1)*2-1,2] = kl
      preLines[(j-1)*2,1] = stoppingPoint
      preLines[(j-1)*2,2] = kl
      nItems = stoppingPoint - startingPoint
      cluster = normalParams[[3]][stoppingPoint]
      clustersKLDivMat[cluster, i] = 
        clustersKLDivMat[cluster, i] + nItems * kl
      nPerCluster[cluster] = nPerCluster[cluster] + nItems
      startingPoint = stoppingPoint + 1
    }
    preLines[(length(changePoints)-1)*2+1,1] = startingPoint
    stoppingPoint = 3000
    kl = calcKLDiv(studentParams[[1]][startingPoint,], studentParams[[2]][startingPoint,], 
                   studentParams[[3]][startingPoint,], estStateProbs, 
                   normalParams[[1]][startingPoint,], normalParams[[2]][startingPoint,], 
                   trueStateProbs, 10000)
    preLines[(length(changePoints)-1)*2+1,2] = kl
    preLines[length(changePoints)*2,1] = stoppingPoint
    preLines[length(changePoints)*2,2] =  kl 
    nItems =  stoppingPoint - startingPoint
    cluster = normalParams[[3]][stoppingPoint]
    clustersKLDivMat[cluster, i] = 
      clustersKLDivMat[cluster, i] + nItems * kl
    nPerCluster[cluster] = nPerCluster[cluster] + nItems
    
    outputLinesFile = paste(outputDir, outputFileString, "KLLines", i, ".csv", sep="")
    write.table(preLines, file=outputLinesFile)
    
    clustersKLDivMat[6, i] = sum(clustersKLDivMat[1:5,i])/sum(nPerCluster)
    for (j in 1:5){
      if(nPerCluster[j] > 1){
        clustersKLDivMat[j,i] = clustersKLDivMat[j,i]/nPerCluster[j]
      }
    }
  }
  return(clustersKLDivMat)
}



trueInputDir = "data/simulatedSequences/"
adaptiveInputDir = "data/resultsFP/"
fittedModelDir = "data/fittedHMMSequences/"
outputDir = "data/toLatexKLDiv/"

set.seed(3)


print("Preprocessing")
ppClustersKLDivMat = matrix(0, nrow = 6, ncol=4)
for (i in 1:4){
  trueStateProbsFile = paste(trueInputDir, "stationaryDistr", i, ".csv", sep="")
  trueStateProbs = read.table(trueStateProbsFile)
  estStateProbsFile = paste(fittedModelDir, "stationaryDistr", i, ".txt", sep="")
  estStateProbs = read.table(estStateProbsFile)
  normalParams = getNormalParams(trueInputDir, i)
  studentParams = getPreprocessingStudentParams(adaptiveInputDir, i)
  changePoints = getAllPreprocessingChangePoints(trueInputDir, adaptiveInputDir, i)
  preLines = matrix(nrow = 2*length(changePoints), ncol= 2)
  nPerCluster = numeric(4)
  for (j in 2:length(changePoints)){
    startingPoint = changePoints[[j-1]]
    stoppingPoint = changePoints[[j]]-1
    preLines[(j-1)*2-1,1] = startingPoint
    kl = calcKLDiv(studentParams[[1]][startingPoint,], studentParams[[2]][startingPoint,], 
                   studentParams[[3]][startingPoint,], estStateProbs, 
                   normalParams[[1]][startingPoint,], normalParams[[2]][startingPoint,], 
                   trueStateProbs, 10000)
    preLines[(j-1)*2-1,2] = kl
    preLines[(j-1)*2,1] = stoppingPoint
    preLines[(j-1)*2,2] = kl
    nItems = stoppingPoint - startingPoint
    cluster = normalParams[[3]][stoppingPoint]
    ppClustersKLDivMat[cluster, i] = 
      ppClustersKLDivMat[cluster, i] + nItems * kl
    nPerCluster[cluster] = nPerCluster[cluster] + nItems
  }
  startingPoint = changePoints[[length(changePoints)]]
  preLines[(length(changePoints)-1)*2+1,1] = startingPoint
  stoppingPoint = 1000
  kl = calcKLDiv(studentParams[[1]][startingPoint,], studentParams[[2]][startingPoint,], 
                 studentParams[[3]][startingPoint,], estStateProbs, 
                 normalParams[[1]][startingPoint,], normalParams[[2]][startingPoint,], 
                 trueStateProbs, 10000)
  preLines[(length(changePoints)-1)*2+1,2] = kl
  preLines[length(changePoints)*2,1] = stoppingPoint
  preLines[length(changePoints)*2,2] =  kl 
  nItems =  stoppingPoint - startingPoint
  cluster = normalParams[[3]][stoppingPoint]
  ppClustersKLDivMat[cluster, i] = 
    ppClustersKLDivMat[cluster, i] + nItems * kl
  nPerCluster[cluster] = nPerCluster[cluster] + nItems

  outputLinesFile = paste(outputDir, "preKLLines", i, ".csv", sep="")
  write.table(preLines, file=outputLinesFile)
  
  
  ppClustersKLDivMat[5, i] = sum(ppClustersKLDivMat[1:4,i])/sum(nPerCluster)
  for (j in 1:4){
    if(nPerCluster[j] > 1){
      ppClustersKLDivMat[j,i] = ppClustersKLDivMat[j,i]/nPerCluster[j]
    }
  }
}
print(ppClustersKLDivMat)


print("Full process")
fpClustersKLDivMat = linesMatrixAdaptive(
  trueInputDir, fittedModelDir, adaptiveInputDir, outputDir, "fullProcess")
print(fpClustersKLDivMat)

adaptiveInputDir = "data/resultsNCM/"

print("No create merge")
ncmClustersKLDivMat = linesMatrixAdaptive(
  trueInputDir, fittedModelDir, adaptiveInputDir, outputDir, "noCreateMerge")
print(ncmClustersKLDivMat)

adaptiveInputDir = "data/resultsSP/"

print("Only switch pre")
spClustersKLDivMat = linesMatrixAdaptive(
  trueInputDir, fittedModelDir, adaptiveInputDir, outputDir, "onlySwitchPre")
print(spClustersKLDivMat)
klDivList = list(ppClustersKLDivMat, fpClustersKLDivMat, ncmClustersKLDivMat, spClustersKLDivMat)
for(i in 1:4){
  outputFormattedKLDivSequence(klDivList, i, outputDir) 
}



  
