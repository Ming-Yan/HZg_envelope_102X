#ifndef H_MAKETREE
#define H_MAKETREE
#include "TTree.h"
TTree* makeTTree_signal(string fname)
{
  // Create ROOT TTree filled with a Gaussian distribution in x and a uniform distribution in y

  Double_t* mH = new Double_t;
  Double_t* weights = new Double_t;
  TTree *T1 = new TTree("T1", "data from ascii file");
  T1->Branch("mH", &mH, "mH/D"); // The branch name is m
  T1->Branch("weights", &weights, "weights/D");
  T1->ReadFile(fname.c_str(), "mH:weights");
  return T1 ;
}
TTree* makeTTree(string fname)
{
  // Create ROOT TTree filled with a Gaussian distribution in x and a uniform distribution in y

  Double_t* mH = new Double_t;
  Double_t* weights = new Double_t;
  TTree *T1 = new TTree("T1", "data from ascii file");
  T1->Branch("mH", &mH, "mH/D"); // The branch name is m
  T1->ReadFile(fname.c_str(), "mH");
  return T1 ;
}
#endif
