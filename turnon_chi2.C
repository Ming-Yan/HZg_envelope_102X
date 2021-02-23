// #include "RooStepBernstein.h"
// #include "RooGaussStepBernstein.h"
#include "RooExponential.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
//#include "boost/program_options.hpp"                                                                             
#include "boost/lexical_cast.hpp"
#include <RooGaussModel.h>
#include <RooTruthModel.h>
#include <RooDecay.h>
#include "RooAddModel.h"
#include <RooNumConvPdf.h>
#include <RooFFTConvPdf.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TFrame.h"
#include "RooFitResult.h"
#include "Maketree.h"
using namespace RooFit;
void turnon_chi2(int cat)
{
  // gSystem->Load("RooStepBernstein_cxx.so");
  // gSystem->Load("RooGaussStepBernstein_cxx.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  
  TTree* DataTree1 = makeTTree(Form("bkg/ptwei_turnon/HiggsMass_ptwei_bkg_ele_mu_cat%d_2020_full.txt",cat));
  int totalev = DataTree1->GetEntriesFast();
  RooRealVar CMS_hzg_mass("CMS_hzg_mass", "CMS_hzg_mass", 100, 170, "GeV") ;
  CMS_hzg_mass.setRange("window",100, 170);
  RooDataSet data(Form("data_obs_ele_mu_cat%d_2020",cat), " ", RooArgSet(CMS_hzg_mass), Import(*DataTree1));
  RooArgList storedPdfs("store");
  // TFile *fout  = new TFile(Form("turnon_cat%d.root",cat),"recreate");
  //generetic Bernstein polynomials
  // RooAbsReal nll;
  // RooMinimizer *bern1_min = new RooMinimizer(nll);
  // bern1_min->minimize("Minuit2","Scan");

  
  RooRealVar mean("mean","mean",0);
  RooRealVar sigma_b1("sigma_b1","sigma_b1",3,3.,10);
  RooRealVar sigma_b2("sigma_b2","sigma_b2",3,3.,10);
  RooRealVar sigma_b3("sigma_b3","sigma_b3",3,3.,10);
  RooRealVar sigma_b4("sigma_b4","sigma_b4",3,3.,10);
  RooRealVar sigma_b5("sigma_b5","sigma_b5",3,3.,10);
  RooRealVar step_b1("step_b1","step_b1",107,105,110);
  RooRealVar step_b2("step_b2","step_b2",107,105,110);
  RooRealVar step_b3("step_b3","step_b3",107,100,110);
  RooRealVar step_b4("step_b4","step_b4",107,100,110);
  RooRealVar step_b5("step_b5","step_b5",107,100,110);
  RooRealVar p0("p0","p0",15);
  RooRealVar b1p1("b1p1","b1p1",0.3,-10.,5);
  RooRealVar b2p1("b2p1","b2p1",0.3,-10.,5);
  RooRealVar b2p2("b2p2","b2p2",3.0819,-5.,5);
  RooRealVar b3p1("b3p1","b3p1",0.3,-5.,5);
  RooRealVar b3p2("b3p2","b3p2",0.3,-5.,5);
  RooRealVar b3p3("b3p3","b3p3",0.3,-5.,5);
  RooRealVar b4p1("b4p1","b4p1",0.3,-5.,5);
  RooRealVar b4p2("b4p2","b4p2",0.3,-5.,5);
  RooRealVar b4p3("b4p3","b4p3",-0.1,-5.,5);//untag
  RooRealVar b4p4("b4p4","b4p4",0.1,-5.,5);//VBF&lepton
  RooRealVar b5p1("b5p1","b5p1",4.5,-5.,5);//untag
  RooRealVar b5p2("b5p2","b5p2",1.67,-5.,5);
  RooRealVar b5p3("b5p3","b5p3",0.24,-5.,5);
  RooRealVar b5p4("b5p4","b5p4",0.18,-5.,5);
  RooRealVar b5p5("b5p5","b5p5",0.17,-5.,5);
  RooGaussStepBernstein bern1("bern1","bern1",CMS_hzg_mass,mean,sigma_b1,step_b1, RooArgList(p0,b1p1));
  RooGaussStepBernstein bern2("bern2","bern2",CMS_hzg_mass,mean,sigma_b2,step_b2, RooArgList(p0,b2p1,b2p2));
  RooGaussStepBernstein bern3("bern3","bern3",CMS_hzg_mass,mean,sigma_b3,step_b3, RooArgList(p0,b3p1,b3p2,b3p3));
  RooGaussStepBernstein bern4("bern4","bern4",CMS_hzg_mass,mean,sigma_b4,step_b4, RooArgList(p0,b4p1,b4p2,b4p3,b4p4));
  RooGaussStepBernstein bern5("bern5","bern5",CMS_hzg_mass,mean,sigma_b5,step_b5, RooArgList(p0,b5p1,b5p2,b5p3,b5p4,b5p5));
  //testing with generic power law   
  //cat1
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",6.7,5.,8.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",108.32,106,109);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.0,-8,-2.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",3.6576e-04,0.0,0.5);
  //cat2
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",4.8,3.,10);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",108,107,110);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-15,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.002 ,0.0,1.);
  //cat3
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",3.5,3.,8.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",106,105,108);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.5,-15,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.02 ,0.0,1.);
  //cat4
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",3.3,2.5,8.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",106.9,105,108);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.3453,-10,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",3.3018e-02 ,0.0,1.);
  //cat501
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",6.,0.1,10);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",108,106,109);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-15,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.0002 ,0.0,1.);
  //cat502
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",6.,0.1,10);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",104.6,104,105);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-8.,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.0006,0,1.);
  //cat503
  RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",4.,0.1,5.);
  RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",104.6,104,105);
  RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.2,-8,-1.);
  RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.001,0.000001,1.);
  //cat6789
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",4.95,1.,8.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",105.7,104,106);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-4.,-8,-1.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.001,0.000001,1.);
  //cat1
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",6.6041,3.,8.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",108.2,107,109);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.0519,-10,-2.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",4.3299e-05,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-4.5113,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",4.0951e-13,0,1.);
  //cat2
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",4.8,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",108,107,110);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.5,-10,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.4,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-5.,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.01,0,1.);
  //cat3
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",3.5,3.,5.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",107.57,106,108);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-8.8928,-10,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.89 ,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-6.5,-10.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.001,0,1.);
  //cat4
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",3.5,3.,8.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",107,105,108);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.3671,-10,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",9.6543e-01,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-4.267,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",1.1001e-06,0,1.);
  //cat501
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",8.505,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",108,107,109);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.4708,-10,-2.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",8.4306e-01,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-6.1246,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",1.4141e-06,0,1.);
  //cat502
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",5.,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",104.6,100,110);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.1,-10.,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.72,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-5.1,-6.,-1.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.0003,0,1.);
  //cat503
  RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",5,3.,10);
  RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",104.6,104,105);
  RooRealVar p1_pow3("p1_pow3","p1_pow3",-5,-8,-3.);
  RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.3,0.,1.);
  RooRealVar p3_pow3("p3_pow3","p3_pow3",-4,-5.,-1.);
  RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.01,0,1.);
  //cat6789
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",4.8,3.,6.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",105.7,104,107);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.1,-10.,-1.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",1.e-06,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-3.0,-6.,-1.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.81,0,1.);
  //cat1
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",6,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108,107,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.04,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.6,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat2
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",4.8,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108,107,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.04,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.6,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat3
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",3.8618,2.,5.);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",107.34,105,108);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-7.0724,-10,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",9.1185e-01,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.4414,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",5.7468e-06,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-8.4821,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",2.8846e-03,0.0,1.);
  //cat4
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",3.5546,2.,5.);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",106.91,106,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-8.8,-11,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",1.e-05,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-6.5,-8.,-3.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",1.9518e-02,0.,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-6.3,-8.,-3.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",2.0522e-02,0.0,1.);
  //cat501
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",6,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108,107,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.04,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.6,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat502
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",6,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",104.6,100,110);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.04,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.6,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat503
  RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",4.8,3.,10);
  RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",104.6,104,105);
  RooRealVar p1_pow5("p1_pow5","p1_pow5",-6,-12,-5);
  RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.07,0.,1.);
  RooRealVar p3_pow5("p3_pow5","p3_pow5",-5.,-10.,-2.);
  RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.03,0,1.);
  RooRealVar p5_pow5("p5_pow5","p5_pow5",-4,-10.,-1.);
  RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.0001,0.001,1.);
  //cat6789
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",4.95,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",105.9,104,107);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-6.2,-12,-4);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.5,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-3.,-10.,0.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.37,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-3.8,-10.,0.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.5,0.0,1.);

  RooGenericPdf step_pow1("step_pow1", "step_pow1", "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(CMS_hzg_mass,turnon_pow1,p1_pow1,cp1_pow1));//step*(ax^b)
  RooGenericPdf step_pow3("step_pow3", "step_pow3", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(CMS_hzg_mass,turnon_pow3,p1_pow3,cp1_pow3,p3_pow3,cp3_pow3));//step*(ax^b+cx^d)
  RooGenericPdf step_pow5("step_pow5", "step_pow5", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(CMS_hzg_mass,turnon_pow5,p1_pow5,cp1_pow5,p3_pow5,cp3_pow5,p5_pow5,cp5_pow5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_pow1("gau_pow1","gau_pow1",CMS_hzg_mass,mean,sigma_pow1);
  RooGaussModel gau_pow3("gau_pow3","gau_pow3",CMS_hzg_mass,mean,sigma_pow3);
  RooGaussModel gau_pow5("gau_pow5","gau_pow5",CMS_hzg_mass,mean,sigma_pow5);
  RooFFTConvPdf gauxpow1("gauxpow1","gauxpow1",CMS_hzg_mass,step_pow1,gau_pow1);
  RooFFTConvPdf gauxpow3("gauxpow3","gauxpow3",CMS_hzg_mass,step_pow3,gau_pow3);
  RooFFTConvPdf gauxpow5("gauxpow5","gauxpow5",CMS_hzg_mass,step_pow5,gau_pow5);
  
  //testing with generic Laurent
  //cat1
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",6,3.,10);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",108,107,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000000001,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99999,0,1.);
  //cat2
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4.8,4.,6.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",108,107,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000001,0.0,0.1);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.95,0,0.5);
  //cat3
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",5.,2.5,8.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",106,103,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",9.5367e-07,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.9999,0,1.);
  //cat4
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",3.1,2.5,8.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",106.5,103,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",4.6222e-10,0.0,0.01);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",4.2351e-01,0.,1.);
  //cat501
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",6,3.,10);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",108,107,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000000001,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99999,0,1.);
  //cat502
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",6,3.,10);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",104.6,100,110);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000000001,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99999,0,1.);
  //cat503
  RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4.9997,3.,8.);
  RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",106.5,104,108);
  RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",8.8602e-01,0.0,1.);//untag
  RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",5.7348e-08,0,1.);
  //cat6789
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4.95,3.,6.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",105.9,103,107);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",2.7933e-01,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",8.3706e-03,0,1.);
  //cat1
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",6,3.,10);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",108,107,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.0000000000000000000000000000000000001,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.9999999999999999999999999999999,0.,1.05);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.0,0,1.);
  //cat2
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.4,4.,5.5);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",108,107,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.2e-10,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.99,0.1,1.);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",1.8e-13,0,1.);
  //cat3
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",3.0238,2.5,5.);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",105.69,105,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.0493e-07,0.,0.1);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",9.8448e-01 ,0.8,1.0);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",9.4361e-08,0,0.1);
  //cat4
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",3.1,1.,8);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",106.6,105,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",3.2949e-06,0.,0.1);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",9.5000e-01,0.9,1.0);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",1.0828e-07,0,0.1);
  //cat501
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",6,3.,10);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",108,107,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.00000000001,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.999999999,0.,1.);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.0,0,1.);
  //cat502
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",6,3.,10);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",104.6,100,110);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.00000000001,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.999999999,0.,1.);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.0,0,1.);
  //cat503
  RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.95,3.,6.);
  RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",105.9,104,107);
  RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.7839e-06,0.,1.);
  RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.95,0.,1.05);
  RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.1,0,1.);
  //cat6789
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.95,3.,6.);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",105.9,104,107);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.7839e-06,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.95,0.,1.05);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.1,0,1.);
 
  //cat1
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",6,3.,10);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",108,107,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.0000000001,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.0000001,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.0,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.999999999,0,1.);
  //cat2
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4.8,3.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",108.,107,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",3.6156e-09,0.,.001);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",4.4677e-05,0.,.1);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",7.5589e-10,0,.0001);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",5.3754e-01,0.3,0.7);
  // cat3
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",3.,1.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",106.5,104,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",1.8872e-15,0.,.001);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",3.2066e-13,0.,.1);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",9.5240e-17,0,.0001);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",9.9935e-01,0.9,1.);
  //cat4
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",3.53,2.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",107.2,105,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",9.3514e-10,0.,.001);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",9.5466e-08,0.,.1);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",8.7939e-11,0,.0001);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",9.0697e-01,0.9,1.);
   //cat501
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",6,3.,10);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",108,107,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.0000000001,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.0000001,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.0,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.999999999,0,1.);
  // cat502
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",6,3.,10);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",104.6,100,110);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.0000000001,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.0000001,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.0,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.999999999,0,1.);
  //cat503
  RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4.,3.,8);
  RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",104.6,104,105);
  RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.05,0.,1.);
  RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.06,0.,1.);
  RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.000001,0,1.);
  RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.02,0,1.);
  //cat6789
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",5.,3.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",106.,104,107);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",3.0862e-02,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",4.1937e-02,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",5.7030e-03,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",7.5681e-01,0,1.);
  RooGenericPdf step_lau1("step_lau1", "step_lau1", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(CMS_hzg_mass,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
  RooGenericPdf step_lau2("step_lau2", "step_lau2", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(CMS_hzg_mass,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
  RooGenericPdf step_lau3("step_lau3", "step_lau3", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(CMS_hzg_mass,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
  RooGaussModel gau_lau1("gau_lau1","gau_lau1",CMS_hzg_mass,mean,sigma_lau1);
  RooGaussModel gau_lau2("gau_lau2","gau_lau2",CMS_hzg_mass,mean,sigma_lau2);
  RooGaussModel gau_lau3("gau_lau3","gau_lau3",CMS_hzg_mass,mean,sigma_lau3);
  RooFFTConvPdf gauxlau1("gauxlau1","gauxlau1",CMS_hzg_mass,step_lau1,gau_lau1);
  RooFFTConvPdf gauxlau2("gauxlau2","gauxlau2",CMS_hzg_mass,step_lau2,gau_lau2);
  RooFFTConvPdf gauxlau3("gauxlau3","gauxlau3",CMS_hzg_mass,step_lau3,gau_lau3);


  // testing with generic exponential
  
  //cat1
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5,3.,10);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 108,107,109);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.03,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat2
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",4.7,3.,6.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",108.,105.,110);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.052,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.9,0,1.);
  //cat3
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",3.85,3.,5.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",107,106.,108.);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.05,-0.5,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.05,0,0.1);
  //cat4
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",3.56,3.,6.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",106.97,105.,108);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-4.8370e-02,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",2.9940e-01,0,1.);
  //cat501
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5,3.,10);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 108,107,109);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.03,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat502
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5,3.,10);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",104.6,100,110);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.04,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat503
  RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",4,3.,10);
  RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 104.6,104,105);
  RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.03,-0.7,0.);
  RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat6789
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5.0,5.,6.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 106.9,105,107.5);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.0228,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.023,0,1.);

  //cat1
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",5.,3.,10);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",108,107,109);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.04,-0.2,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.3,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.03,-0.2,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.1,0,1.);
  //cat2
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",4.8,3.,10);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",108.,100,110);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.07,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.75,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.04,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.5,0,1.);
  //cat3
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",3.8507,3.,5.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",106.95,105,108.);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-5.4503e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",9.6120e-01,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-4.8149e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",1.8388e-07,0,.1);
  //cat4
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",3.58,3.,5.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",106.7,105,108);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-3.4577e-02 ,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",1.2078e-02,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-6.5367e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",5.7985e-01,0,1.);
  //cat501
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",8.104,6.,10);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",108,107,109);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-4.6408e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",2.3801e-05,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-4.6408e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",2.3053e-04,0,1.);
  //cat502
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",5.1,3.,10);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",104.6,100,110);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.006,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.75,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.5,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.6,0,1.);
  //cat503
  RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",3.832,2.,6.);
  RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",104.6,104,105);
  RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.025,-0.1,0.);
  RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.17,0.,1.);
  RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.024,-0.1,0.);
  RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",1.1120e-06,0,1.);
  //cat6789
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",4.95,3.,8.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",105.9,104,107);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-3.4369e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.42,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-1.8018e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.10,0,1.);
  
  //cat1
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.5,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",108,107,109);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.15,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.09 ,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.1,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.05,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.1,0,1.);
  //cat2
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",4.8,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",108,100,110);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.175,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.05,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.2,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.00001,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.000001,0,1.);
  //cat3
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",3.5,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",106,100,110);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.175,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.05,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.2,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.00001,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.000001,0,1.);
  //cat4
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",3.5,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",106,100,110);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.175,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.05,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.2,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.00001,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.000001,0,1.);
  //cat501
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.5,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",108,107,109);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-1.8351e-01,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",1.5511e-08,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-1.0083e-01 ,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",7.9829e-04,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-4.6365e-02,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",1.6977e-01,0,1.);
  //cat502
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.5,3.,10);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",104.6,100,110);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.054,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.0004,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.042 ,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.06,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.035,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.02,0,1.);
  //cat503
  RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.,3.,10);
  RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",104.6,104,105);
  RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.2,-0.5,0.);
  RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.03,-0.5,0.);
  RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.2,0,1.);
  RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.00001,-0.5,0.);
  RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.000001,0,1.);
  //cat6789
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.,3.,7.);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",106,104,106.5);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.035,-0.6,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.4,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.02,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.1,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.01,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.002,0,1.);
  RooGenericPdf step_exp1("step_exp1", "step_exp1", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(CMS_hzg_mass,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
  RooGenericPdf step_exp3("step_exp3", "step_exp3", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(CMS_hzg_mass,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
  RooGenericPdf step_exp5("step_exp5", "step_exp5", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(CMS_hzg_mass,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_exp1("gau_exp1","gau_exp1",CMS_hzg_mass,mean,sigma_exp1);
  RooGaussModel gau_exp3("gau_exp3","gau_exp3",CMS_hzg_mass,mean,sigma_exp3);
  RooGaussModel gau_exp5("gau_exp5","gau_exp5",CMS_hzg_mass,mean,sigma_exp5);
  RooFFTConvPdf gauxexp1("gauxexp1","gauxexp1",CMS_hzg_mass,step_exp1,gau_exp1);
  RooFFTConvPdf gauxexp3("gauxexp3","gauxexp3",CMS_hzg_mass,step_exp3,gau_exp3);
  RooFFTConvPdf gauxexp5("gauxexp5","gauxexp5",CMS_hzg_mass,step_exp5,gau_exp5);
  
  
	int stat=1;
  int ntries=0;
	double minnll=10e8;

  RooArgSet *params_pow1 = gauxpow1.getParameters((const RooArgSet*)(0));
  RooArgSet *params_pow3 = gauxpow3.getParameters((const RooArgSet*)(0));
  RooArgSet *params_pow5 = gauxpow5.getParameters((const RooArgSet*)(0));
  RooArgSet *params_exp1 = gauxexp1.getParameters((const RooArgSet*)(0));
  RooArgSet *params_exp3 = gauxexp3.getParameters((const RooArgSet*)(0));
  RooArgSet *params_exp5 = gauxexp5.getParameters((const RooArgSet*)(0));
  RooArgSet *params_lau1 = gauxlau1.getParameters((const RooArgSet*)(0));
  RooArgSet *params_lau2 = gauxlau2.getParameters((const RooArgSet*)(0));
  RooArgSet *params_lau3 = gauxlau3.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern1 = bern1.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern2 = bern2.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern3 = bern3.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern4 = bern4.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern5 = bern5.getParameters((const RooArgSet*)(0));
  RooFitResult *pow1_fit;  RooFitResult *exp1_fit;  RooFitResult *lau1_fit;
  RooFitResult *pow3_fit;  RooFitResult *exp3_fit;  RooFitResult *lau2_fit;
  RooFitResult *pow5_fit;  RooFitResult *exp5_fit;  RooFitResult *lau3_fit;
  RooFitResult *bern1_fit; RooFitResult *bern2_fit; RooFitResult *bern3_fit; RooFitResult *bern4_fit; RooFitResult *bern5_fit;
	// while (stat!=0){
	//   if (ntries>=5) break;
	  
  //   stat = pow1_fit->status();
	//   minnll = pow1_fit->minNll();
	//   if (stat!=0) params_pow1->assignValueOnly(pow1_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat!=0){
	// //   if (ntries>=5 ) break;
    // pow1_fit = gauxpow1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","Scan"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // pow3_fit = gauxpow3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // pow5_fit = gauxpow5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // exp1_fit = gauxexp1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // exp3_fit = gauxexp3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // exp5_fit = gauxexp5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // lau1_fit = gauxlau1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // lau2_fit = gauxlau2.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // lau3_fit = gauxlau3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // // // // // bern1_fit = bern1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  bern2_fit = bern2.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(1)); //FIXME
	  // bern3_fit = bern3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // bern4_fit = bern4.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // bern5_fit = bern5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  

  storedPdfs.add(gauxexp1);
  storedPdfs.add(gauxexp3);
  storedPdfs.add(gauxexp5);
  storedPdfs.add(gauxpow1);
  storedPdfs.add(gauxpow3);
  storedPdfs.add(gauxpow5);
  storedPdfs.add(gauxlau1);
  storedPdfs.add(gauxlau2);
  storedPdfs.add(gauxlau3);
  storedPdfs.add(bern1);
  storedPdfs.add(bern2);
  storedPdfs.add(bern3);
  storedPdfs.add(bern4);
  storedPdfs.add(bern5);
  
  // RooWorkspace *ws =  new RooWorkspace();
  // RooCategory catIndex("catIndex","c");
  // RooMultiPdf *pdf = new RooMultiPdf("CMS_hzg_bkgshape","all pdfs",catIndex,storedPdfs);
  // RooRealVar nBackground("CMS_hzg_bkgshape_norm","nbkg",data.sumEntries(),0,10E8);
  // ws->SetName("multipdf");
  // ws->import(*pdf);
  // ws->import(nBackground);
  // ws->import(catIndex);
  // ws->import(data);
  // fout->cd();
  // ws->Write();
  // fout->Close();
  RooPlot* xframe1  = CMS_hzg_mass.frame() ;
  CMS_hzg_mass.setRange("blind1",100,120) ;
  CMS_hzg_mass.setRange("blind2",130,170);
  data.plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(80),CutRange("blind2")) ;
  gauxpow1.plotOn(xframe1,RooFit::Name("gauxpow1"),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashed));
  gauxpow3.plotOn(xframe1,RooFit::Name("gauxpow3"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashed));
  gauxpow5.plotOn(xframe1,RooFit::Name("gauxpow5"),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashed));
  gauxexp1.plotOn(xframe1,RooFit::Name("gauxexp1"),LineColor(TColor::GetColor("#FF8589")),LineStyle(kDotted));
  gauxexp3.plotOn(xframe1,RooFit::Name("gauxexp3"),LineColor(TColor::GetColor("#81C4EF")),LineStyle(kDotted));
  gauxexp5.plotOn(xframe1,RooFit::Name("gauxexp5"),LineColor(TColor::GetColor("#FFD35C")),LineStyle(kDotted));
  gauxlau1.plotOn(xframe1,RooFit::Name("gauxlau1"),LineColor(TColor::GetColor("#A30005")),LineStyle(kDashDotted));
  gauxlau2.plotOn(xframe1,RooFit::Name("gauxlau2"),LineColor(TColor::GetColor("#136090")),LineStyle(kDashDotted));
  gauxlau3.plotOn(xframe1,RooFit::Name("gauxlau3"),LineColor(TColor::GetColor("#E0A500")),LineStyle(kDashDotted));
  // bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(TColor::GetColor("#FF595E")));
  bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(TColor::GetColor("#FFCA3A")));
  bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(TColor::GetColor("#8AC926")));
  bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(TColor::GetColor("#1982C4")));
  bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(TColor::GetColor("#6A4C93")));

    xframe1->SetMinimum(0.0001);
  xframe1->Draw();
  TLegend* leg4 = new TLegend(0.6,0.6,0.9,0.9);
  leg4->SetNColumns(3);
  leg4->AddEntry(xframe1->findObject("gauxpow1"), "pow1", "l");
  leg4->AddEntry(xframe1->findObject("gauxpow3"), "pow3", "l");
  leg4->AddEntry(xframe1->findObject("gauxpow5"), "pow5", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp1"), "exp1", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp3"), "exp3", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp5"), "exp5", "l");
  leg4->AddEntry(xframe1->findObject("gauxlau1"), "lau1", "l");
  leg4->AddEntry(xframe1->findObject("gauxlau2"), "lau2", "l");
  leg4->AddEntry(xframe1->findObject("gauxlau3"), "lau3", "l");
  // leg4->AddEntry(xframe1->findObject("bern1"), "bern1", "l");
  leg4->AddEntry(xframe1->findObject("bern2"), "bern2", "l");
  leg4->AddEntry(xframe1->findObject("bern3"), "bern3", "l");
  leg4->AddEntry(xframe1->findObject("bern4"), "bern4", "l");
  leg4->AddEntry(xframe1->findObject("bern5"), "bern5", "l");
    leg4->AddEntry(xframe1->findObject("data"), "data", "lep");

  leg4->Draw("same");
  gPad->Print(Form("cat%d_turn.pdf",cat));
  cout<<"=======>Category"<<cat<<endl;
  
//  cout<<"Bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();
//   if(bern2_fit->status()!=0){cout<<"Bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();}
//   if(bern3_fit->status()!=0){cout<<"Bern3:"<<bern3_fit->status()<<endl;bern3_fit->Print();}
//   if(bern4_fit->status()!=0){cout<<"Bern4:"<<bern4_fit->status()<<endl;bern4_fit->Print();}
//   if(bern5_fit->status()!=0){cout<<"Bern5:"<<bern5_fit->status()<<endl;bern5_fit->Print();}
//   if(pow1_fit->status()!=0){cout<<"Pow1:"<<pow1_fit->status()<<endl;pow1_fit->Print();}
//   if(pow3_fit->status()!=0){cout<<"Pow3:"<<pow3_fit->status()<<endl;pow3_fit->Print();}
//   if(pow5_fit->status()!=0){cout<<"Pow5:"<<pow5_fit->status()<<endl;pow5_fit->Print();}
//  if(lau1_fit->status()!=0){cout<<"Lau1:"<<lau1_fit->status()<<endl;lau1_fit->Print();}
//   if(lau2_fit->status()!=0){cout<<"Lau2:"<<lau2_fit->status()<<endl;lau2_fit->Print();}
//   if(lau3_fit->status()!=0){cout<<"Lau3:"<<lau3_fit->status()<<endl;lau3_fit->Print();}
//   if(exp1_fit->status()!=0){cout<<"Exp1:"<<exp1_fit->status()<<endl;exp1_fit->Print();}
//   if(exp3_fit->status()!=0){cout<<"Exp3:"<<exp3_fit->status()<<endl;exp3_fit->Print();}
//   if(exp5_fit->status()!=0){cout<<"Exp5:"<<exp5_fit->status()<<endl;exp5_fit->Print();}
}
