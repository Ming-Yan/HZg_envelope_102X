#!/bin/bash
ulimit -s unlimited
set -e
cd /afs/cern.ch/user/m/milee/CMSSW_10_2_13/src/flashggFinalFit/Plots/SplusBModels1/toys
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700 
eval `scramv1 runtime -sh`

#itoy=0
for itoy in $(seq 0 99)
do
#Generate command
combine /afs/cern.ch/user/m/milee/CMSSW_10_2_13/src/flashggFinalFit/Plots/cat1_post_snap.root -m 125.380 -M GenerateOnly --saveWorkspace --toysFrequentist --bypassFrequentistFit -t 1 --setParameters r=0.214 -s -1 -n _${itoy}_gen_step --snapshotName MultiDimFit

#Fit command
mv higgsCombine_${itoy}_gen_step*.root gen_${itoy}.root
root -l -b -q getsnap.C\(\"gen\"\,\"1\"\,${itoy}\)
combine gen_${itoy}.root -m 125.380 -M MultiDimFit -P r --floatOtherPOIs=1 --saveWorkspace --toysFrequentist --bypassFrequentistFit -t 1 --setParameters r=0.214 -s -1 -n _${itoy}_fit_step --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --snapshotName "MultiDimFit"

#Throw command
mv higgsCombine_${itoy}_fit_step*.root fit_${itoy}.root
root -l -b -q getsnap.C\(\"fit\"\,\"1\"\,${itoy}\)
combine fit_${itoy}.root -m 125.380 --snapshotName MultiDimFit -M GenerateOnly --saveToys --toysFrequentist --bypassFrequentistFit -t -1 -n _${itoy}_throw_step --setParameters r=0

mv higgsCombine_${itoy}_throw_step*.root toy_${itoy}.root
rm gen_${itoy}.root fit_${itoy}.root
done
