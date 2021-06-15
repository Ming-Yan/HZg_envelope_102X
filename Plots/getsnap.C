void getsnap(string toy,string cat="1",int i=0){
  TFile *fin = TFile::Open(Form("/afs/cern.ch/user/m/milee/CMSSW_10_2_13/src/flashggFinalFit/Plots/cat%s_post_snap.root",cat.c_str()));
  //for (int  i = 0; i < 1;i++){
    TFile *f = TFile::Open(Form("%s_%d.root",toy.c_str(),i),"update");
  RooWorkspace *w = (RooWorkspace*)fin->Get("w");
  RooArgSet *snap = (RooArgSet*)f->Get("MultiDimFit");
  w->saveSnapshot("MultiDimFit", *snap,true);
  w->expensiveObjectCache().clearAll(); w->expensiveObjectCache().clearAll(); w->writeToFile(Form("%s_%d.root",toy.c_str(),i));
  //}
}

