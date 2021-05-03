
export("calcSuffStats")
calcSuffStats <- function(execTimesOccProbs, start, stop, nStates){
  A0 <- rep(0, nStates)
  A1 <- rep(0, nStates)
  A2 <- rep(0, nStates)
  n <- stop - start
  for (i in 1:nStates){
    #    print(dataFrame$viterbiStateProbs[, i+1])
    A0[i] <- A0[i] + sum(execTimesOccProbs[start:stop, i+1])
    A1[i] <- A1[i] + sum(execTimesOccProbs[start:stop, i+1]*execTimesOccProbs[start:stop, nStates+2])
    A2[i] <- A2[i] + sum(execTimesOccProbs[start:stop, i+1]*execTimesOccProbs[start:stop, nStates+2]^2)
  }
  return(list(A0, A1, A2, n))
}






