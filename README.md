# FLASHgg Final Fits for H->Zg analysis 
The Final Fits package is a series of scripts which are used to run the final stages of the CMS Hgg analysis: signal modelling, background modelling, datacard creation and final statistical interpretation and final result plots.

## Pdf model adpated with turnon fit

[x] Bernstein 

[x] Power law

[x] Laurent

[ ] Exponential

-> Power law and Laurent series are done by discrete convolution

-> the exponential can use analytic way


## Download and setup instructions

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git cms-init

# Install the GBRLikelihood package which contains the RooDoubleCBFast implementation
git clone git@github.com:jonathon-langford/HiggsAnalysis.git
# Install Combine as per the documentation here: cms-analysis.github.io/HiggsAnalysis-CombinedLimit/
git clone git@github.com:cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

# Compile external libraries
cd HiggsAnalysis
cmsenv
scram b -j 9

# Install Flashgg Final Fit packages
cd ..
git clone -b dev_runII_102x git@github.com:cms-analysis/flashggFinalFit.git
cd flashggFinalFit/
```

Two packages need to be built with their own makefiles, if needed. Please note that there will be verbose warnings from BOOST etc, which can be ignored. So long as the `make` commands finish without error, then the compilation happened fine.:

```
cd ${CMSSW_BASE}/src/flashggFinalFit/Background
make
```

## Contents
The FLASHgg Finals Fits package contains several subfolders which are used for the following steps:


* Create the Background Model (see `Background` dir)
* Generate a Datacard (see `Datacard` dir)
* Run `combine` and generate statistical interpretation plots. (see `Plots/FinalResults` dir)

Each of the relevant folders are documented with specific `README.md` files.

## Known issues

Recently some issues with memory have been observed with the workspaces (probably because there are so many processes and tags now). Crashes can occur due to a `std::bad_alloc()` error, which for now I have managed to circumvent by submitting to the batch (this is at Imperial College), e.g. for making the photon systematic dat files and the S+B fits. The problem is due to all the workspaces being loaded by the WSTFileWrapper class, so at some point this should be revisited and improved somwhow. 

## Updates in dev_runII_102x branch

* New, easier to navigate submission scripts: option for config file
* Integration with HTCondor
* Python module for replacement dataset map: used when too few entries to construct signal model
* Pruning: removes processes below threshold in datacard



The modes are used for the following (run in sequential order):
  * `datacard` - build the .txt datacard using the S & B models. The yield variations from systematics are also calculated and specified in the datacard. To merge datacards for different years then use the `combineCards.py` script (in combine).
  * `combine`  - compile the RooWorkspace from the .txt datacard. Run the fit in combine. Input options are specified in `Plots/FinalResults/combineHarvesterOptions_Template.dat`
  * `combinePlots` - create plots from finished combine jobs. Options are specified in `Plots/FinalResults/combinePlotsOptions_Template.dat`

