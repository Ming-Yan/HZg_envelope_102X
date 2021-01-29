for c in VBF  #nm1_all untag VBF     #2 3 4 501 502 503 6789 #nm1_all untag VBF
do
#    python makeToys.py --inputWSFile higgsCombineacceff_${c}_bf.MultiDimFit.mH125.38.root --loadSnapshot MultiDimFit --ext ${c} --queue workday  --batch condor
 #   cd  SplusBModels${c}/toys/jobs/
   # condor_submit sub_toys.sub
  #  cd /afs/cern.ch/user/m/milee/CMSSW_10_2_13/src/flashggFinalFit/Plots
     python makeSplusBModelPlot_new.py --inputWSFile higgsCombineacceff_${c}_bf.MultiDimFit.mH125.38.root --loadSnapshot MultiDimFit --cats all  --doBkgRenormalization 1 --doToyVeto 1 --ext ${c} --nBins 55
done
