# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 15:13:31 2020

@author: afe02
"""

import numpy as np
from numpy.random import seed
import csv
from BOSplitTraceHelperBayesianForwardBackward import LikelihoodsCalculator
from BOSplitTraceHelperBayesianForwardBackward import recursiveBO
import matplotlib.pyplot as plt


seed(12345)

filesPrefix = "../data/fittedHMMSequences/"
outputFilesPrefix = "../data/splitIndicesSequencesPreprocess/"

GLR_limits = [-20, -20, -20, -20]
for i in range(1,5):
    print("sequence " + str(i))
    timesProbsFile = filesPrefix + "executionTimesOccupancyProbabilities" + str(i) + ".csv"   
    statProbsFile = filesPrefix + "stationaryDistr" +  str(i) + ".txt"
    normalParamsFile = filesPrefix + "normalParams" + str(i) + ".txt"
    transitionMatrixFile = filesPrefix + "transitionMatrix" + str(i) + ".txt"
    data = np.genfromtxt(timesProbsFile, delimiter=',', skip_header=1)
    weights = np.genfromtxt(statProbsFile, delimiter=' ', skip_header=0)
    normalParams = np.genfromtxt(normalParamsFile, delimiter=' ', skip_header=0)
    transitionMatrix = np.genfromtxt(transitionMatrixFile, delimiter=' ', skip_header=0)
    nStates = data.shape[1] - 2
    length = 999
    likelihoodCalc = LikelihoodsCalculator(0, length, data[:, 1:nStates+1], \
                                           data[:, nStates+1], nStates, length, 5, \
                                           normalParams, 20*weights, transitionMatrix, weights)
    # can be removed
    n = 80
    x = np.random.rand(n).reshape(n,1)
    logLikelihoodSums = likelihoodCalc.negLogLikelihoodSumVec(x)
    plt.plot(length*x, logLikelihoodSums, 'x')

    x = range(0, len(data[:, nStates+1]))
    plt.figure()
    plt.plot(x, data[:, nStates+1], 'x')
    # to here - kept to keep identical results with random methods

    result = recursiveBO(0, length, 20, likelihoodCalc, length, GLR_limits[i-1])

    # can be removed    
    print(result)
    plt.figure()
    plt.plot(x, data[:, nStates+1], 'x')
    for r in result:
        plt.axvline(x=length*r)
    # to here - kept to keep identical results with random methods
    
    splitIndicesFileName = outputFilesPrefix + "splitIndicesSeq" + str(i) + ".csv"
    with open(splitIndicesFileName, 'w', newline='') as f:
        indicesWriter = csv.writer(f, delimiter=',')
        indicesWriter.writerow(("index",))
        for j in range(len(result)):
            indicesWriter.writerow((int(result[j][0]*length),))


