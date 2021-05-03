# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 15:13:31 2020

@author: afe02
"""

import GPy
import GPyOpt
from scipy.special import loggamma
from scipy.stats import t
import numpy as np
import math

def initializeHyperParametersNormalGamma(normalParams, nPseudoObs, nStates):
    hyperParams = np.zeros((nStates, 4))
    hyperParams[:, 0] = normalParams[:, 0]
    hyperParams[:, 1] = nPseudoObs
    hyperParams[:, 2] = nPseudoObs/2
    hyperParams[:, 3] = hyperParams[:, 2]*normalParams[:,1]*normalParams[:,1]
    return hyperParams

def updateHyperparamsNormalGamma(a0, a1, a2, hyperParams, nStates, threshold = 0.001):
  validIndices = np.where(a0 > threshold)
  hyperParams_out = np.zeros((nStates, 4))
  hyperParams_out[:,:] = hyperParams[:,:]
  mu_in = hyperParams[:,0]
  nu_in = hyperParams[:,1]
  # alpha is redundant - nu/2
  alpha_in = hyperParams[:,2]
  beta_in = hyperParams[:,3]
  # mu update
  hyperParams_out[validIndices,0] = (nu_in[validIndices]*mu_in[validIndices] + a1[validIndices])/(nu_in[validIndices] + a0[validIndices])
  # nu update
  hyperParams_out[validIndices, 1] = nu_in[validIndices] + a0[validIndices]
  # alpha is redundant - nu/2
  hyperParams_out[validIndices,2] = alpha_in[validIndices] + a0[validIndices]/2
  # beta update
  mu_diff = a1[validIndices]/a0[validIndices] - mu_in[validIndices]
  hyperParams_out[validIndices, 3] = beta_in[validIndices] + \
    0.5*(a2[validIndices] - a1[validIndices]*a1[validIndices]/a0[validIndices] + \
           nu_in[validIndices] * a0[validIndices] * \
           mu_diff*mu_diff/(nu_in[validIndices] + a0[validIndices]))
  return hyperParams_out

def studentParams(hyperParams, nStates):
  studentParams = np.zeros((nStates, 3))
  means = hyperParams[:,0]
  scaling = np.sqrt(hyperParams[:,3]*(hyperParams[:,1] + 1)/(hyperParams[:,1]*hyperParams[:,2]))
  dof = hyperParams[:,1]
  studentParams[:,0] = means
  studentParams[:,1] = scaling
  studentParams[:,2] = dof
  return studentParams


def forwardBackwardStudent(executionTimes, nStates, studentParams, initProbs, transitionMatrix):
  alphas = np.zeros((len(executionTimes), nStates))
  scaling = np.zeros((len(executionTimes), 1))
  betas = np.zeros((len(executionTimes), nStates))
  gammas = np.zeros((len(executionTimes), nStates))
  studentStateProbs = np.zeros((len(executionTimes), nStates))
  for i in range(len(executionTimes)):
      for j in range(nStates):
          studentStateProbs[i,j] = t.pdf(executionTimes[i], studentParams[j, 2], \
                           loc=studentParams[j, 0], scale=studentParams[j, 1])
  for i in range(nStates):
    alphas[0, i] = initProbs[i]*studentStateProbs[0,i]
    betas[len(executionTimes)-1,i] = 1
  scaling[0, 0] = 1/np.sum(alphas[0, :])
  alphas[0,] = alphas[0,:]*scaling[0,0]
  logLike = math.log(scaling[0,0])
  for i in range(1, len(executionTimes)):
    for j in range(nStates):
      probj = np.sum(transitionMatrix[:,j]*alphas[i-1,])
      alphas[i,j] = probj*studentStateProbs[i,j]
    scaling[i,0] =  1/np.sum(alphas[i,:])
    alphas[i,:] = alphas[i,:]*scaling[i,0]
    logLike -= math.log(scaling[i,0])
  gammas[len(executionTimes)-1,:] = alphas[len(executionTimes)-1,:]
  for i in range(len(executionTimes)-2, 0, -1):
      for j in range(nStates):
          betas[i,j] = np.sum(transitionMatrix[j,:]*betas[i+1,:]*studentStateProbs[i+1,:])
      betas[i,:] = betas[i,:]*scaling[i,0]
      gammas[i,:] = betas[i,:]*alphas[i,:]
      gammas[i,:] = gammas[i,:]/np.sum(gammas[i,:])
  return (logLike, scaling, gammas)


    

def logLikelihoodPart(hyperParams):
#    mu = hyperParams[:,0]
    nu = hyperParams[:,1]
    alpha = hyperParams[:,2]
    beta = hyperParams[:,3]
    res = - loggamma(alpha) + alpha*np.log(beta) + 0.5*np.log(nu)
    return res


def suffStatsHelper(start, x, stop, summationArray, nStates):
    intFrac = np.modf(x)
    suffStats = np.zeros((2, nStates))
    suffStats[0, :] = summationArray[start:int(intFrac[1])+1, :].sum(axis=0) 
    suffStats[0, :] += intFrac[0]*summationArray[int(intFrac[1])+1, :]
    suffStats[1, :] = summationArray[int(intFrac[1])+1:stop, :].sum(axis=0)
    suffStats[1, :] += (1-intFrac[0])*summationArray[int(intFrac[1]), :]
    return suffStats

# todo clean a0limit
class LikelihoodsCalculator:
    def __init__(self, start, stop, occupancyProbabilities, executionTimes, nStates, \
                 xScale, a0Limit, normalParams, nPseudoObs, transitionMatrix, initProbs):
        self.mStart = start
        self.mStop = stop
        self.mOccupancyProbabilities = occupancyProbabilities
        self.mA1Matrix = (occupancyProbabilities.T * executionTimes).T
        self.mA2Matrix = (occupancyProbabilities.T * (executionTimes*executionTimes)).T
        self.mNStates = nStates
        self.mXScale = xScale
        self.mA0Limit = a0Limit
        self.mExecutionTimes = executionTimes
        self.mTransitionMatrix = transitionMatrix
        self.mInitProbs = initProbs
        self.mInitHyperParams = initializeHyperParametersNormalGamma(normalParams, nPseudoObs, nStates)
        
    
          
    def a0s(self, x):
        return suffStatsHelper(self.mStart, x, self.mStop, self.mOccupancyProbabilities, self.mNStates)
   
    def a1s(self, x):
        return suffStatsHelper(self.mStart, x, self.mStop, self.mA1Matrix, self.mNStates)

    def a2s(self, x):
        return suffStatsHelper(self.mStart, x, self.mStop, self.mA2Matrix, self.mNStates)
    
    def negLogLikelihoodSum(self, x):
        A0s = self.a0s(x)
        A1s = self.a1s(x)
        A2s = self.a2s(x)
        hyperParamsLow = updateHyperparamsNormalGamma(A0s[0,:], A1s[0, :], A2s[0, :], self.mInitHyperParams, self.mNStates)
        hyperParamsHigh = updateHyperparamsNormalGamma(A0s[1,:], A1s[1, :], A2s[1, :], self.mInitHyperParams, self.mNStates)
        studentLow = studentParams(hyperParamsLow, self.mNStates)
        studentHigh = studentParams(hyperParamsHigh, self.mNStates)
        etLow = self.mExecutionTimes[self.mStart:int(x)]
        etHigh = self.mExecutionTimes[int(x)+1:self.mStop]
        lowFBStudent = forwardBackwardStudent(etLow, self.mNStates, studentLow, self.mInitProbs, self.mTransitionMatrix)
        highFBStudent = forwardBackwardStudent(etHigh, self.mNStates, studentHigh, self.mInitProbs, self.mTransitionMatrix)
        return -(lowFBStudent[0] + highFBStudent[0])
    
    def negLogLikelihoodSumVec(self, xVec):
        negLogLikelihoodSums = np.zeros(len(xVec))
        for i in range(0, len(xVec)):
            negLogLikelihoodSums[i] = self.negLogLikelihoodSum(self.mXScale*xVec[i])
        return negLogLikelihoodSums
    
    def GLR(self, x):
        a0sAll = self.mOccupancyProbabilities[self.mStart:self.mStop, :].sum(axis=0)
        a1sAll = self.mA1Matrix[self.mStart:self.mStop, :].sum(axis=0)
        a2sAll = self.mA2Matrix[self.mStart:self.mStop, :].sum(axis=0)
        hyperParamsAll = updateHyperparamsNormalGamma(a0sAll, a1sAll, a2sAll, self.mInitHyperParams, self.mNStates)
        studentAll = studentParams(hyperParamsAll, self.mNStates)
        etAll = self.mExecutionTimes[self.mStart:self.mStop]
        FBStudentAll = forwardBackwardStudent(etAll, self.mNStates, studentAll, self.mInitProbs, self.mTransitionMatrix)
        logLikelihoodAll = FBStudentAll[0]
        negllsx = self.negLogLikelihoodSum(x*self.mXScale)
        return logLikelihoodAll + negllsx


def recursiveBO(start, stop, minLength, likelihoodCalc, scaleX, GLR_limit):
    if (stop - start < minLength*2):
        return []
    likelihoodCalc.mStart = start
    likelihoodCalc.mStop = stop
    bounds_likelihood = [{'name': 'var_1', 'type': 'continuous', 'domain': ((start + minLength)/scaleX, (stop - minLength)/scaleX)}] 
    k_likelihood = GPy.kern.RBF(1)
    BO_g = GPyOpt.methods.BayesianOptimization(f = likelihoodCalc.negLogLikelihoodSumVec,           
                                               domain = bounds_likelihood,       
                                               kernel = k_likelihood,          
                                               acquisition='EI',
                                               acquisition_par = 0.1)
    BO_g.model.noise_var = 0.001
    BO_g.run_optimization(max_iter=40)

    split_x = BO_g.x_opt
    GLR_x = likelihoodCalc.GLR(split_x)
    if (GLR_x > GLR_limit):
        return []
    print(split_x*scaleX)
    print(GLR_x)
    split_left = recursiveBO(start, int(split_x*scaleX), minLength, likelihoodCalc, scaleX, GLR_limit)
    split_right = recursiveBO(int(split_x*scaleX), stop, minLength, likelihoodCalc, scaleX, GLR_limit)
    return split_left + [split_x] + split_right
    
    
    
