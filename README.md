# AdaptiveETBayes
Code related to Adaptive Runtime Estimate of Task Execution Times using Bayesian Modeling

# R-related version info
RStudio Version 1.3.1093

R version 4.0.3 (2020-10-10)

platform       x86_64-w64-mingw32          
arch           x86_64                      
os             mingw32                     
system         x86_64, mingw32             
svn rev        79318                       
nickname       Bunny-Wunnies Freak Out    
packages:
depmixS4 version 1.4-2
data.tree version 1.0.0
MASS version 7.3.53
ggplot version 3.3.2

# Python-related info
Python and GpyOpt used for preprocessing step 2 (find points of model change)

Install anaconda.

In anaconda prompt:
```console
>conda update -y anaconda 
>conda update -y numpy scipy matplotlib
>conda update -y jupyter
```

Install environment python_3_gpyopt with dependencies:
```console
conda env create -f environment.yml
```

# Generate simulated sequences
Source from RStudio in AdaptiveETBayes project

createEvalSequences/simulateMarkovClusterSequences.R

Result found in data/simulatedSequences

# Preprocessing step 1: Fit HMM
Run from RStudio in AdaptiveETBayes project

preprocessing1fitHMM/simMarkovFitHMMSequences.R

Result found in data/fittedHMMSequences

# Preprocessing step 2: Find points of model change
In anaconda prompt:
```console
conda activate python_3_gpyopt
cd preprocessing2pointsModelChange
python BOSplitTraceSimulatedBayesianSequences.py
```

Result found in data/splitIndicesSequencesPreprocess

# Preprocessing step 3: Cluster segments + adaptive step
Adaptive step, full process (FP)

Source from RStudio in AdaptiveETBayes project

preprocessing3clusterAdaptive/clusterSequencesGLRSlidingWindowFP.R

Result found in data/resultsFP

Adaptive step, no create/ merge (NCM)

Source from RStudio in AdaptiveETBayes project

preprocessing3clusterAdaptive/clusterSequencesGLRSlidingWindowNCM.R

Result found in data/resultsNCM

Adaptive step, switch preprocessing (SP)

Source from RStudio in AdaptiveETBayes project

preprocessing3clusterAdaptive/clusterSequencesGLRSlidingWindowSP.R

Result found in data/resultsSP

# Generate sequential lines for ground truth and estimated distributions (Fig 3)
Source from RStudio in AdaptiveETBayes project

evaluationFigures/exportDistrLines.R

Result found in 

data/toLatexDistrLinesFP

data/toLatexDistrLinesNCM

data/toLatexDistrLinesSP

# Calculate KL divergence and generate lines for graphs (Fig 3)
Source from RStudio in AdaptiveETBayes project

evaluationFigures/KLDivCalcSequence.R

Result found in

data/toLatexKLDiv


