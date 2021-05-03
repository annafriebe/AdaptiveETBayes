


outputTrueLines <- function(trueInputDir, indexIn, outputDir, indexOut){
  trueSplittingPointsFile = paste(trueInputDir, "splittingPointsClustersSequence", indexIn, ".csv", sep="")
  trueSplittingPoints = read.table(trueSplittingPointsFile)
  trueMeansFile = paste(trueInputDir, "meansSequence", indexIn, ".csv", sep="")
  trueMeans = read.table(trueMeansFile)
  trueStddevsFile = paste(trueInputDir, "stdDevsSequence", indexIn, ".csv", sep="")
  trueStdDevs = read.table(trueStddevsFile)
  lines = matrix(nrow = 2*nrow(trueSplittingPoints), ncol=10)
  for (row in 1:(nrow(trueSplittingPoints)-1)){
    startIndex = (row-1)*2 + 1
    lines[startIndex,1] = trueSplittingPoints[row,1]
    lines[startIndex,2] = trueMeans[trueSplittingPoints[row,2],1]
    lines[startIndex,3] = trueMeans[trueSplittingPoints[row,2],2]
    lines[startIndex,4] = trueMeans[trueSplittingPoints[row,2],3]
    lines[startIndex,5] = trueMeans[trueSplittingPoints[row,2],1] + 
      trueStdDevs[trueSplittingPoints[row,2],1]
    lines[startIndex,6] = trueMeans[trueSplittingPoints[row,2],2] + 
      trueStdDevs[trueSplittingPoints[row,2],2]
    lines[startIndex,7] = trueMeans[trueSplittingPoints[row,2],3] + 
      trueStdDevs[trueSplittingPoints[row,2],3]
    lines[startIndex,8] = trueMeans[trueSplittingPoints[row,2],1] - 
      trueStdDevs[trueSplittingPoints[row,2],1]
    lines[startIndex,9] = trueMeans[trueSplittingPoints[row,2],2] - 
      trueStdDevs[trueSplittingPoints[row,2],2]
    lines[startIndex,10] = trueMeans[trueSplittingPoints[row,2],3] - 
      trueStdDevs[trueSplittingPoints[row,2],3]
    lines[startIndex+1,1] = trueSplittingPoints[row+1,1]
    lines[startIndex+1,2:10] = lines[startIndex,2:10]
  }
  startIndex = (nrow(trueSplittingPoints)-1)*2 + 1
  lines[startIndex, 1] = trueSplittingPoints[nrow(trueSplittingPoints),1]
  lines[startIndex,2] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],1]
  lines[startIndex,3] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],2]
  lines[startIndex,4] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],3]
  lines[startIndex,5] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],1]+ 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],1]
  lines[startIndex,6] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],2]+ 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],2]
  lines[startIndex,7] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],3]+ 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],3]
  lines[startIndex,8] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],1]- 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],1]
  lines[startIndex,9] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],2]- 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],2]
  lines[startIndex,10] = trueMeans[trueSplittingPoints[nrow(trueSplittingPoints),2],3]- 
    trueStdDevs[trueSplittingPoints[nrow(trueSplittingPoints), 2],3]
  lines[startIndex + 1,1] = 3000
  lines[startIndex+1,2:10] = lines[startIndex,2:10]
  outputLinesFile = paste(outputDir, "lines", indexOut, ".csv", sep="")
  write.table(lines, file=outputLinesFile)
  
}

outputPreProcessingLines <- function(preprocessingInputDir, indexIn, outputDir, indexOut){
  changingPointsClusters = paste(preprocessingInputDir, "clusterIndexSequence", indexIn, ".csv", sep="")
  pointsClusters = read.table(changingPointsClusters)
  nClusters = max(pointsClusters[,2])
  studentParams = list()
  for (cluster in 1:nClusters){
    studentParamsFileName = paste(preprocessingInputDir, "studentParams", indexIn, "_", cluster, ".txt", sep="")
    studentParams[[cluster]] = read.table(studentParamsFileName)
  }
  nStates = nrow(studentParams[[1]])
  preLines = matrix(nrow = 2*nrow(pointsClusters), ncol= 1+3*nStates)
  for(row in 1:(nrow(pointsClusters)-1)){
    startIndex = (row-1)*2 + 1
    preLines[startIndex, 1] = pointsClusters[row,1]
    studentParamsCluster = studentParams[[pointsClusters[row,2]]]
    for (state in 1:nStates){
      mean = studentParamsCluster[state, 1]
      scale = studentParamsCluster[state, 2]
      dof = studentParamsCluster[state, 3]
      stdDev = scale*sqrt(dof/(dof-2))
      preLines[startIndex, (state-1)*3 + 2] = mean
      preLines[startIndex, (state-1)*3 + 3]= mean + stdDev
      preLines[startIndex, (state-1)*3 + 4]= mean - stdDev
    }
    preLines[startIndex+1,1] = pointsClusters[row+1, 1]
    preLines[startIndex+1, 2:(1+3*nStates)] = 
      preLines[startIndex, 2:(1+3*nStates)]
  }
  startIndex = (nrow(pointsClusters)-1)*2 + 1
  preLines[startIndex, 1] = pointsClusters[nrow(pointsClusters), 1]
  studentParamsCluster = studentParams[[pointsClusters[nrow(pointsClusters),2]]]
  for (state in 1:nStates){
    mean = studentParamsCluster[state, 1]
    scale = studentParamsCluster[state, 2]
    dof = studentParamsCluster[state, 3]
    stdDev = scale*sqrt(dof/(dof-2))
    preLines[startIndex, (state-1)*3 + 2] = mean
    preLines[startIndex, (state-1)*3 + 3]= mean + stdDev
    preLines[startIndex, (state-1)*3 + 4]= mean - stdDev
  }
  preLines[startIndex+1,1] = 1000
  preLines[startIndex+1, 2:(1+3*nStates)] = 
    preLines[startIndex, 2:(1+3*nStates)]
  
  outputLinesFile = paste(outputDir, "preLines", indexOut, ".csv", sep="")
  write.table(preLines, file=outputLinesFile)
}

outputAdaptiveLines <- function(adaptiveInputDir, indexIn, outputDir, indexOut){
  adaptiveStudentChangesFile = paste(adaptiveInputDir, "adaptiveStudentParamsChangeIndices", indexIn, ".csv", sep="")
  adaptiveChangePoints = read.table(adaptiveStudentChangesFile, header=FALSE, stringsAsFactors= FALSE, sep=",")
  adaptiveStudentParamsFile = paste(adaptiveInputDir, "adaptiveStudentParams", indexIn, ".csv", sep="")
  adaptiveStudentParams = read.table(adaptiveStudentParamsFile, header=FALSE, stringsAsFactors= FALSE, sep=",")
  nStates = nrow(adaptiveStudentParams)/(3*nrow(adaptiveChangePoints))

  adaptiveLines = matrix(nrow = 2*nrow(adaptiveChangePoints), ncol= 1+3*nStates)
  for(row in 1:(nrow(adaptiveChangePoints)-1)){
    startIndex = (row-1)*2 + 1
    adaptiveLines[startIndex, 1] = adaptiveChangePoints[row,1]
    startReadIndexParams = (row-1)*(3*nStates)
    meansIndex = startReadIndexParams
    for (state in 1:nStates){
      meansIndex = meansIndex + 1
      scaleIndex = meansIndex + nStates
      dofIndex = scaleIndex + nStates
      means = adaptiveStudentParams[meansIndex,1]
      dof = adaptiveStudentParams[dofIndex,1]
      stdDev = adaptiveStudentParams[scaleIndex,1]*sqrt(dof/(dof-2))
      adaptiveLines[startIndex, (state-1)*3 + 2] = means
      adaptiveLines[startIndex, (state-1)*3 + 3]= means + stdDev
      adaptiveLines[startIndex, (state-1)*3 + 4]= means - stdDev
    }
    adaptiveLines[startIndex+1,1] = adaptiveChangePoints[row+1, 1]
    adaptiveLines[startIndex+1, 2:(1+3*nStates)] = 
      adaptiveLines[startIndex, 2:(1+3*nStates)]
  }
  startIndex = (nrow(adaptiveChangePoints)-1)*2 + 1
  adaptiveLines[startIndex, 1] = adaptiveChangePoints[nrow(adaptiveChangePoints), 1]
  startReadIndexParams = (nrow(adaptiveChangePoints)-1)*3*nStates
  meansIndex = startReadIndexParams
  for (state in 1:nStates){
    meansIndex = meansIndex + 1
    scaleIndex = meansIndex + nStates
    dofIndex = scaleIndex + nStates
    
    means = adaptiveStudentParams[meansIndex,1]
    dof = adaptiveStudentParams[dofIndex,1]
    stdDev = adaptiveStudentParams[scaleIndex,1]*sqrt(dof/(dof-2))
    adaptiveLines[startIndex, (state-1)*3 + 2] = means
    adaptiveLines[startIndex, (state-1)*3 + 3]= means + stdDev
    adaptiveLines[startIndex, (state-1)*3 + 4]= means - stdDev
  }
  
  adaptiveLines[startIndex+1,1] = 3000
  adaptiveLines[startIndex+1, 2:(1+3*nStates)] = 
    adaptiveLines[startIndex, 2:(1+3*nStates)]

  outputLinesFile = paste(outputDir, "adaptiveLines", indexOut, ".csv", sep="")
  write.table(adaptiveLines, file=outputLinesFile)
  
}

trueInputDir = "data/simulatedSequences/"
adaptiveInputDir = "data/resultsFP/"
outputDir = "data/toLatexDistrLinesFP/"


for (i in 1:4){
  outputTrueLines(trueInputDir, i, outputDir, i)
  outputPreProcessingLines(adaptiveInputDir, i, outputDir, i)
  outputAdaptiveLines(adaptiveInputDir, i, outputDir, i)
}

trueInputDir = "data/simulatedSequences/"
adaptiveInputDir = "data/resultsNCM/"
outputDir = "data/toLatexDistrLinesNCM/"


for (i in 1:4){
  outputTrueLines(trueInputDir, i, outputDir, i+10)
  outputPreProcessingLines(adaptiveInputDir, i, outputDir, i+10)
  outputAdaptiveLines(adaptiveInputDir, i, outputDir, i+10)
}

trueInputDir = "data/simulatedSequences/"
adaptiveInputDir = "data/resultsSP/"
outputDir = "data/toLatexDistrLinesSP/"


for (i in 1:4){
  outputTrueLines(trueInputDir, i, outputDir, i+20)
  outputPreProcessingLines(adaptiveInputDir, i, outputDir, i+20)
  outputAdaptiveLines(adaptiveInputDir, i, outputDir, i+20)
}

