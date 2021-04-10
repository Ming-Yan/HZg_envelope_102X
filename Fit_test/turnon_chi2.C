
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
// #include "ModGaus11.h"
using namespace RooFit;
void turnon_chi2(int cat)
{
  // gSystem->Load("RooStepBernstein_cxx.so");
  // gSystem->Load("RooGaussStepBernstein_cxx.so");
  //gSystem->Load("/afs/cern.ch/work/m/milee/MYcode/limit/ModGaus11_cxx.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  // gSystem->Load("ModGaus11_cxx.so");
  TTree* DataTree1 = makeTTree(Form("bkg/ptwei_turnon/HiggsMass_ptwei_bkg_ele_mu_cat%d_2020_full.txt",cat));
  int totalev = DataTree1->GetEntriesFast();
  RooRealVar CMS_hzg_mass("CMS_hzg_mass", "CMS_hzg_mass", 100, 170, "GeV") ;
  CMS_hzg_mass.setRange("window",100, 170);
  RooDataSet data(Form("data_obs_ele_mu_cat%d_2020",cat), " ", RooArgSet(CMS_hzg_mass), Import(*DataTree1));
  RooDataHist datahist(Form("datahist_ele_mu_cat%d_2020",cat),Form("datahist_ele_mu_cat%d_2020",cat),CMS_hzg_mass,data);
  RooArgList storedPdfs("store");
  TFile *fout  = new TFile(Form("newturnon_cat%d.root",cat),"recreate");


  //generetic Bernstein polynomials
  RooRealVar mean("mean","mean",0);
  RooRealVar sigma_b1("sigma_b1","sigma_b1",4.,3.,10);
  RooRealVar sigma_b2("sigma_b2","sigma_b2",4.,3.,10);
  RooRealVar sigma_b3("sigma_b3","sigma_b3",4.,3.,10);
  RooRealVar sigma_b4("sigma_b4","sigma_b4",4.,3.,10);
  RooRealVar sigma_b5("sigma_b5","sigma_b5",4.,3.,8.);
  RooRealVar step_b1("step_b1","step_b1",107,100,110);
  RooRealVar step_b2("step_b2","step_b2",107,100,110);
  RooRealVar step_b3("step_b3","step_b3",107,100,110);
  RooRealVar step_b4("step_b4","step_b4",105,100,110);
  RooRealVar step_b5("step_b5","step_b5",105,103,108);
  RooRealVar p0("p0","p0",15);
  RooRealVar b1p1("b1p1","b1p1",0.3,-25.,25);
  RooRealVar b2p1("b2p1","b2p1",0.3,-25.,25);
  RooRealVar b2p2("b2p2","b2p2",3.,-25.,25);
  RooRealVar b3p1("b3p1","b3p1",0.3,-25.,25);
  RooRealVar b3p2("b3p2","b3p2",0.3,-25.,25);
  RooRealVar b3p3("b3p3","b3p3",0.3,-25.,25.);
  RooRealVar b4p1("b4p1","b4p1",0.3,-25.,25.);
  RooRealVar b4p2("b4p2","b4p2",0.3,-25.,25.);
  RooRealVar b4p3("b4p3","b4p3",-0.1,-25.,25.);//untag
  RooRealVar b4p4("b4p4","b4p4",0.1,-25.,25);//VBF&lepton
  RooRealVar b5p1("b5p1","b5p1",4.5,-15.,15.);//untag
  RooRealVar b5p2("b5p2","b5p2",1.67,-15.,15.);
  RooRealVar b5p3("b5p3","b5p3",0.24,-15.,15.);
  RooRealVar b5p4("b5p4","b5p4",0.18,-15.,15.);
  RooRealVar b5p5("b5p5","b5p5",0.17,-15.,15.);
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
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",108.5,107,110);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-15,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.002 ,0.0,1.);
  //cat3
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",4,3.,8.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",106.6,106,107);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.5,-15,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.02 ,0.0,1.);
  //cat4
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",3.5717,3.,4.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",106.9,106,108);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.3453,-7,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",4.7131e-05 ,0.0,1.);
  //cat501
  RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",6.,0.1,10);
  RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",108,106,109);
  RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-15,-5.);
  RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.0002 ,0.0,1.);
  //cat502
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",6.,0.1,10);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",104.6,104,105);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.6,-8.,-5.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.0006,0,1.);
  //cat503
  // RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",4.,0.1,5.);
  // RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",104.6,104,105);
  // RooRealVar p1_pow1("p1_pow1","p1_pow1",-5.2,-8,-1.);
  // RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.001,0.000001,1.);
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
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",5.0,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",108.5,107,110);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-7.,-10,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",7.7276e-01,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-6,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",6.3834e-05,0,1.);
  //cat3
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",4.,3.,5.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",106.6,106,107);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3", -7.0864,-9,-6.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.99,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-4.8364,-6.,-3.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",7.1789e-14,0,0.001);
  //cat4
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",3.5,3.,8.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",107,105,108);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.3671,-10,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",9.6543e-01,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-4.267,-8.,-2.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",1.1001e-06,0,1.);
  //cat501
  RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",8.5830,7.,9.);
  RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",108.47,108,109);
  RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.5,-10,-5.);
  RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",5.3658e-03 ,0.,1.);
  RooRealVar p3_pow3("p3_pow3","p3_pow3",-6.,-8.,-5.);
  RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",4.4629e-10,0,0.001);
  //cat502
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",5.,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",104.6,100,110);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.1,-10.,-5.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.72,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-5.1,-6.,-1.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.0003,0,1.);
  //cat503
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",4.,3.,10);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",104.6,104,105);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-5,-8,-3.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.3,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-4,-5.,-1.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.01,0,1.);
  //cat6789
  // RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",4.8,3.,6.);
  // RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",105.7,104,107);
  // RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.1,-10.,-1.);
  // RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",1.e-06,0.,1.);
  // RooRealVar p3_pow3("p3_pow3","p3_pow3",-3.0,-6.,-1.);
  // RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.81,0,1.);
  //cat1
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",6.6,5.,7.);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108.42,107,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-5.2664,-7,-4);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5", 2.6960e-06,0.,.01);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-9.7909,-12.,-8.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",3.3933e-04,0,.1);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-6.0477,-8.,-5.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",5.5285e-02,0.,0.1);
  //cat2
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",5.,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108.64,107,109);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-10.4,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",7.4192e-01,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.0937,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",9.0955e-01,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-6.2482,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",2.1590e-03,0.001,1.);
  //cat3
//   RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",4.,3.,5.);
//   RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",106.6,106,107);
//  RooRealVar p1_pow5("p1_pow5","p1_pow5",-8.8,-11,-5);
//   RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",1.e-05,0.,1.);
//   RooRealVar p3_pow5("p3_pow5","p3_pow5",-6.5,-8.,-3.);
//   RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",1.9518e-02,0.,1.);
//   RooRealVar p5_pow5("p5_pow5","p5_pow5",-6.3,-8.,-3.);
//   RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",2.0522e-02,0.0,1.);
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
  RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",6,3.,10);
  RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",108,107,109);
  RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.04,0.,1.);
  RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.6,-10.,-1.);
  RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat502
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",5.,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",104.6,100,110);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-6,-15,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.1,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-7.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.2,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.001,0.001,1.);
  //cat503
  // RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",4.,3.,10);
  // RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",104.6,104,105);
  // RooRealVar p1_pow5("p1_pow5","p1_pow5",-6,-12,-5);
  // RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.07,0.,1.);
  // RooRealVar p3_pow5("p3_pow5","p3_pow5",-5.,-10.,-2.);
  // RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.03,0,1.);
  // RooRealVar p5_pow5("p5_pow5","p5_pow5",-4,-10.,-1.);
  // RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.0001,0.00001,1.);
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
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4.05,4.,6.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",108,107,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",5.8468e-08,0.0,0.1);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.25,0,1.);
  //cat3
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4.,3.,6.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",107,106.,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",9.5367e-07,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.9999,0,1.);
  //cat4
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",3.1,2.5,8.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",106.5,103,109);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",4.6222e-10,0.0,0.01);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",4.2351e-01,0.,1.);
  //cat501
  RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",6,3.,10);
  RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",108,107,109);
  RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000000001,0.0,1.);//untag
  RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99999,0,1.);
  //cat502
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",6,3.,10);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",104.6,100,110);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.000000001,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99999,0,1.);
  //cat503
  // RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",3.9,3.,5.);
  // RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",104.8,104,106);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",9.3937e-01,0.0,1.);//untag
  // RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",4.9829e-06,0,1.);
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
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.0,4.,5.5);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",108,107,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.2503e-06,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",9.7549e-01,0.1,1.);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",2.1416e-09,0,1.);
  //cat3
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.,3.,5.);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",107.,106.,108.);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.0493e-07,0.,0.1);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",9.8448e-01 ,0.5,1.0);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",9.4361e-08,0,0.1);
  //cat4
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",3.1,1.,8);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",106.6,105,109);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",3.2949e-06,0.,0.1);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",9.5000e-01,0.9,1.0);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",1.0828e-07,0,0.1);
  //cat501
  RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",6,3.,10);
  RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",108,107,109);
  RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.00000000001,0.,1.);
  RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.999999999,0.,1.);
  RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.0,0,1.);
  //cat502
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",6,3.,10);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",104.6,100,110);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.00000000001,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.999999999,0.,1.);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.0,0,1.);
  //cat503
  // RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4.,3.,6.);
  // RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",104.6,104,106);
  // RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",1.7839e-06,0.,1.);
  // RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.95,0.,1.05);
  // RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.1,0,1.);
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
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",107.78,107,109);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",5.e-09,0.,.001);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",4.e-08,0.,.1);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",7e-11,0,.0001);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",5.1e-01,0.,1.);
  // cat3
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4.,1.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",107,104,109);
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
  RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",6,3.,10);
  RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",108,107,109);
  RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.0000000001,0.,1.);
  RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.0000001,0.,1.);
  RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.0,0,1.);
  RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.999999999,0,1.);
  // cat502
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4.8,3.,6);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",105,100,110);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",9.5716e-11,0.,0.000001);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.0000001,0.,0.01);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.0,0,0.000001);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",1.,0,1.);
  //cat503
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4.,3.,8);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",104.6,104,105);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.05,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.06,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.000001,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.02,0,1.);
  //cat6789
  // RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",5.,3.,6.);
  // RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",106.,104,107);
  // RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",3.0862e-02,0.,1.);
  // RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",4.1937e-02,0.,1.);
  // RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",5.7030e-03,0,1.);
  // RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",7.5681e-01,0,1.);
  
  //cat1
  // RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",6,3.,10);
  // RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",108,107,109);
  // RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",0.0001,0.,0.01);
  // RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",0.01,0.,0.1);
  // RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",3.2487e-08,0,0.000001);
  // RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",6.6948e-01,0.5,1.);
  // RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",2.3252e-11 ,0,0.1);
  //cat2
//   RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",4.8,3.,6.);
//   RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",107.78,107,109);
//   RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",2.1580e-07,0.,.0001);
//   RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",9.5006e-04,0.,.01);
//   RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",7e-11,0,.0001);
//   RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",6.3471e-01,0.,1.);
//  RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",1.9395e-12,0,.00001);

  // cat3
//   RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",4.,1.,6.);
//   RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",107,104,109);
//   RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",6.6334e-09,0.,.001);
// RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",1.6890e-03,0.,.1);
// RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",7e-11,0,.0001);
// RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",5.9884e-01,0.,1.);
// RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",7.7440e-12,0,0.000001);
  //cat4
  // RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",3.53,2.,6.);
  // RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",107.2,105,109);
  // RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",9.3514e-05,0.,.001);
  // RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",9.5466e-02,0.,.1);
  // RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",8.7939e-17,0,.0001);
  // RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",0.99,0.9,1.);
  // RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",0.,0,1.);
   //cat501
  RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",6,3.,10);
  RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",108,107,109);
  RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",0.0000000001,0.,.0001);
  RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",0.0000001,0.,0.01);
  RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",0.0,0,0.0000001);
  RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",0.999999999,0,1.);
  RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",0.,0,.000001);
  // cat502
//   RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",5.,3.,5.5);
//   RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",106.1,105,108);
//  RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",1.4302e-07,0.,.001);
//   RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",9.5957e-03,0.,.1);
//   RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",7.9832e-10,0,.000001);
//   RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",9.9291e-01,0.99,1.);
//   RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",0.,0,.000000000001);
  //cat503
  // RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",4.,3.,8);
  // RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",104.6,104,105);
  // RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",0.05,0.,1.);
  // RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",0.06,0.,1.);
  // RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",0.000001,0,1.);
  // RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",0.02,0,1.);
  //   RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",0.,0,1.);
  //cat6789
  // RooRealVar sigma_lau4("sigma_lau4","sigma_lau4",5.,3.,6.);
  // RooRealVar turnon_lau4("turnon_lau4","turnon_lau4",106.,104,107);
  // RooRealVar cl1_lau4("cl1_lau4","cl1_lau4",3.0862e-02,0.,1.);
  // RooRealVar cl2_lau4("cl2_lau4","cl2_lau4",4.1937e-02,0.,1.);
  // RooRealVar cl3_lau4("cl3_lau4","cl3_lau4",5.7030e-03,0,1.);
  // RooRealVar cl4_lau4("cl4_lau4","cl4_lau4",7.5681e-01,0,1.);
  // RooRealVar cl5_lau4("cl5_lau4","cl5_lau4",0.,0,1.);
  //cat1
  // RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",6,3.,7);
  // RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",108,107,109);
  // RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",0.0000000001,0.,0.0000001);
  // RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",0.0000001,0.,0.0001);
  // RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",0.0,0,0.00000001);
  // RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",0.5,0,1.);
  // RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",3.3896e-12,0,.1);
  // RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.99,0,1.);
    //cat2
//   RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",4.8,3.,6.);
//   RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",107.78,107,109);
//   RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",5.e-09,0.,.001);
//   RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",4.e-08,0.,.1);
//   RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",7e-11,0,.0001);
//   RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",5.1e-01,0.,1.);
// RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",0.00000,0,0.00001);
// RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.00000,0,0.00001);
  // cat3
//   RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",5.,3.,6.);
//   RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",107,104,109);
//   RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",5.e-09,0.,.001);
// RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",4.e-08,0.,.1);
// RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",7e-11,0,.0001);
// RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",5.1e-01,0.,1.);
// RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",0.0000001,0,1.);
// RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.0000001,0,1.);
  //cat4
  // RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",3.4,2.,6.);
  // RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",106.5,105,109);
  // RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",4.3088e-11,0.,.001);
  // RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",9.7697e-02,0.,.5);
  // RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",4.4971e-12,0,.0001);
  // RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",9.9980e-01 ,0.5,1.);
  // RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",3.6526e-14,0,1.);
  // RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",2.5235e-02,0.0,0.5);
   //cat501
  RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",6,3.,10);
  RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",108,107,109);
  RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",0.0000000001,0.,0.001);
  RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",0.0000001,0.,0.01);
  RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",0.0,0,0.01);
  RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",0.4,0,1.);
  RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",0.,0,0.0000001);
  RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.999999999,0,1.);
  // cat502
  // RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",5,3.,10);
  // RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",106.1,105,108);
  // RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",9.5897e-06,0.,0.001);
  // RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",9.9863e-03,0.,0.1);
  // RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",9.9971e-07,0,0.0001);
  // RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",1.0000e-10,0,0.00000001);
  // RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",0.01,0,.1);
  // RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.999999999,0,1.);
  //cat503
  // RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",4.,3.,8);
  // RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",104.6,104,105);
  // RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",0.05,0.,1.);
  // RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",0.06,0.,1.);
  // RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",0.000001,0,1.);
  // RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",0.02,0,1.);
  //   RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",0.,0,1.);
  // RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",0.999999999,0,1.);
  //cat6789
  // RooRealVar sigma_lau5("sigma_lau5","sigma_lau5",5.,3.,6.);
  // RooRealVar turnon_lau5("turnon_lau5","turnon_lau5",106.,104,107);
  // RooRealVar cl1_lau5("cl1_lau5","cl1_lau5",7.1069e-01 ,0.,1.);
  // RooRealVar cl2_lau5("cl2_lau5","cl2_lau5",4.1937e-02,0.,1.);
  // RooRealVar cl3_lau5("cl3_lau5","cl3_lau5",1.3526e-01,0,1.);
  // RooRealVar cl4_lau5("cl4_lau5","cl4_lau5",6.8827e-04,0,.1);
  // RooRealVar cl5_lau5("cl5_lau5","cl5_lau5",9.4479e-08,0,.1);
  // RooRealVar cl6_lau5("cl6_lau5","cl6_lau5",9.8845e-02,0,1.);

 RooGenericPdf step_lau1("step_lau1", "step_lau1", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(CMS_hzg_mass,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
  RooGenericPdf step_lau2("step_lau2", "step_lau2", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(CMS_hzg_mass,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
  RooGenericPdf step_lau3("step_lau3", "step_lau3", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(CMS_hzg_mass,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
  RooGenericPdf step_lau4("step_lau4", "step_lau4", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2))", RooArgList(CMS_hzg_mass,turnon_lau4,cl1_lau4,cl2_lau4,cl3_lau4,cl4_lau4,cl5_lau4));//step*(ax^b+cx^d)
  RooGenericPdf step_lau5("step_lau5", "step_lau5", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2)+@7*(@0)^(-7))", RooArgList(CMS_hzg_mass,turnon_lau5,cl1_lau5,cl2_lau5,cl3_lau5,cl4_lau5,cl5_lau5,cl6_lau5));//step*(ax^b+cx^d)

  RooGaussModel gau_lau1("gau_lau1","gau_lau1",CMS_hzg_mass,mean,sigma_lau1);
  RooGaussModel gau_lau2("gau_lau2","gau_lau2",CMS_hzg_mass,mean,sigma_lau2);
  RooGaussModel gau_lau3("gau_lau3","gau_lau3",CMS_hzg_mass,mean,sigma_lau3);
  RooGaussModel gau_lau4("gau_lau4","gau_lau4",CMS_hzg_mass,mean,sigma_lau4);
RooGaussModel gau_lau5("gau_lau5","gau_lau5",CMS_hzg_mass,mean,sigma_lau5);

  RooFFTConvPdf gauxlau1("gauxlau1","gauxlau1",CMS_hzg_mass,step_lau1,gau_lau1);
  RooFFTConvPdf gauxlau2("gauxlau2","gauxlau2",CMS_hzg_mass,step_lau2,gau_lau2);
  RooFFTConvPdf gauxlau3("gauxlau3","gauxlau3",CMS_hzg_mass,step_lau3,gau_lau3);
RooFFTConvPdf gauxlau4("gauxlau4","gauxlau4",CMS_hzg_mass,step_lau4,gau_lau4);
RooFFTConvPdf gauxlau5("gauxlau5","gauxlau5",CMS_hzg_mass,step_lau5,gau_lau5);


  // testing with generic exponential
  //cat1
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",6.5,5.,8);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 107.5,107,109);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-4.3265e-02,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",5.1163e-01,0,1.);
  //cat2
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5.0,3.,6.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",108.,105.,110);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.052596,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.98,0,1.);
  //cat3
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",4.,3.,5.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",107,106.,108.);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.05,-0.5,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.05,0,0.1);
  //cat4
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",3.56,3.,4.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",106.97,106.,108);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-4e-02,-0.1,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.01,0,1.);
  //cat501
  RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5,3.,10);
  RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 108,107,109);
  RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.03,-0.7,0.);
  RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat502
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5,3.,10);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",104.6,100,110);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.04,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat503
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",4,3.,10);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 104.6,104,105);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.03,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  //cat6789
  // RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5.0,5.,6.);
  // RooRealVar turnon_exp1("turnon_exp1","turnon_exp1", 106.3,105,106.5);
  // RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.0228,-0.7,0.);
  // RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.023,0,1.);

  //cat1
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",6.2293,6.,7);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",107.83,107,109);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-5.0820e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",1.3882e-02,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3", -4.2712e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",7.2016e-02,0,.5);
  //cat2
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",5.02,3.,6.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",108.41,100,110);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-6.6444e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",6.7480e-01,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-4.1990e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",3.1965e-02,0,.5);
  //cat3
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",4.,3.,5.5);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",107,105,108.);
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
  RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",8.104,6.,10);
  RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",108,107,109);
  RooRealVar p1_exp3("p1_exp3","p1_exp3",-4.6408e-02,-0.5,0.);
  RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",2.3801e-05,0.,1.);
  RooRealVar p3_exp3("p3_exp3","p3_exp3",-4.6408e-02,-0.5,0.);
  RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",2.3053e-04,0,1.);
  //cat502
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",5.1,3.,6);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",105.6,104.6,106);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3", -0.04,-0.1,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",6.3160e-03,0.000001,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-4.2208e-02 ,-1.,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3", 1.6215e-02,0.0000001,1.);
  //cat503
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",3.832,2.,6.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",104.6,104,105);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.025,-0.1,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.17,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.024,-0.1,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",1.1120e-06,0,1.);
  //cat6789
  // RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",4.95,3.,6.);
  // RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",106.3,104,107);
  // RooRealVar p1_exp3("p1_exp3","p1_exp3",-4.8536e-02,-0.5,0.);
  // RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.43,0.,1.);
  // RooRealVar p3_exp3("p3_exp3","p3_exp3",-1.8018e-02,-0.5,0.);
  // RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",3.7871e-02,0,1.);
  
  //cat1
//   RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",6.22,6,7);
//   RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",107.8,107,108);
//  RooRealVar p1_exp5("p1_exp5","p1_exp5",-5.0820e-02,-1.,0.);
//   RooRealVar cp1_exp5("cp1_exp5","cp1_exp5", 4.6641e-02,0.,1.);
//   RooRealVar p3_exp5("p3_exp5","p3_exp5",-4.2712e-02,-0.5,0.);
//   RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",3.7630e-02,0,0.1);
//   RooRealVar p5_exp5("p5_exp5","p5_exp5",-6.8897e-08,-0.1,0.);
//   RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",1.0000e-05,0,.01);
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",8.1943,7.,9);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",107.6,107,109);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-4.6807e-02,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",1.2757e-02 ,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-1.3939e-01 ,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",2.5680e-07,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-4.6554e-02,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",2.4790e-01,0,1.);
  //cat2
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.0658,4.,6.);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",108.3,107,109);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-5.4155e-02,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",6.7451e-01,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-5.5821e-02,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",4.3393e-07,0,.001);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-2.4270e-02,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",2.8593e-04,0,.1);
  //cat3
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",3.85,3.,6);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",107,105,110);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.054,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",8.5202e-01,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.055,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",2.9097e-01,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.059,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",1.9429e-09,0,1.);
  //cat4
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",3.52,3.,6);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",106.8,106,108);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-1.1158e-01,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",9.9456e-01,0.5,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-4.9674e-02,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",1.6118e-02,0,.1);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-5.1149e-05,-0.01,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",1.1229e-06,0,.1);
  //cat501
  RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",8.1943,7.,9);
  RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",107.6,107,109);
  RooRealVar p1_exp5("p1_exp5","p1_exp5",-4.6807e-02,-0.5,0.);
  RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",1.2757e-02 ,0.,1.);
  RooRealVar p3_exp5("p3_exp5","p3_exp5",-1.3939e-01 ,-0.5,0.);
  RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",2.5680e-07,0,1.);
  RooRealVar p5_exp5("p5_exp5","p5_exp5",-4.6554e-02,-0.5,0.);
  RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",2.4790e-01,0,1.);
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
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",3.8366,3.,5.);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",104.6,104,105);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-1.5800e-01,-0.5,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",5.9251e-04,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-2.4977e-02,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",1.6553e-03,0,1.);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-5.7189e-03,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",6.6356e-07 ,0,1.);
  //cat6789
  // RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",4.56,3.,6.);
  // RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",106.3,104,107);
  // RooRealVar p1_exp5("p1_exp5","p1_exp5",-2.7150e-02,-0.6,0.);
  // RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",4.4506e-01,0.,1.);
  // RooRealVar p3_exp5("p3_exp5","p3_exp5",-2.4988e-02,-0.5,0.);
  // RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",3.7690e-06,0,.1);
  // RooRealVar p5_exp5("p5_exp5","p5_exp5",-7.9924e-06,-0.5,0.);
  // RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",2.1386e-03,0,.1);
  RooGenericPdf step_exp1("step_exp1", "step_exp1", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(CMS_hzg_mass,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
  RooGenericPdf step_exp3("step_exp3", "step_exp3", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(CMS_hzg_mass,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
  RooGenericPdf step_exp5("step_exp5", "step_exp5", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(CMS_hzg_mass,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_exp1("gau_exp1","gau_exp1",CMS_hzg_mass,mean,sigma_exp1);
  RooGaussModel gau_exp3("gau_exp3","gau_exp3",CMS_hzg_mass,mean,sigma_exp3);
  RooGaussModel gau_exp5("gau_exp5","gau_exp5",CMS_hzg_mass,mean,sigma_exp5);
  RooFFTConvPdf gauxexp1("gauxexp1","gauxexp1",CMS_hzg_mass,step_exp1,gau_exp1);
  RooFFTConvPdf gauxexp3("gauxexp3","gauxexp3",CMS_hzg_mass,step_exp3,gau_exp3);
  RooFFTConvPdf gauxexp5("gauxexp5","gauxexp5",CMS_hzg_mass,step_exp5,gau_exp5);
  
  //Generic Gaussian
  // RooRealVar turnon_ggau("turnon_ggau","turnon_ggau",112,108,120);
  // RooRealVar turnon_gv0("turnon_gv0","turnon_gv0",1.6,0.,5.);
  // RooRealVar turnon_gv1("turnon_gv1","turnon_gv1",1e-10,0.,0.5);
  // RooRealVar turnon_gs0("turnon_gs0","turnon_gs0",1.5,1.,10);
  // RooRealVar turnon_gs1("turnon_gs1","turnon_gs1",28,5.,50.);
  // ModGaus11 ggau("ggau","ggau",CMS_hzg_mass,turnon_ggau,turnon_gv0,turnon_gv1,turnon_gs0,turnon_gs1);
  
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
  RooArgSet *params_lau4 = gauxlau4.getParameters((const RooArgSet*)(0));
  RooArgSet *params_lau5 = gauxlau5.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern1 = bern1.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern2 = bern2.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern3 = bern3.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern4 = bern4.getParameters((const RooArgSet*)(0));
  RooArgSet *params_bern5 = bern5.getParameters((const RooArgSet*)(0));
  // RooArgSet *params_ggau = ggau.getParameters((const RooArgSet*)(0));
  RooFitResult *pow1_fit;  RooFitResult *exp1_fit;  RooFitResult *lau1_fit;
  RooFitResult *pow3_fit;  RooFitResult *exp3_fit;  RooFitResult *lau2_fit;
  RooFitResult *pow5_fit;  RooFitResult *exp5_fit;  RooFitResult *lau3_fit;
  RooFitResult *lau4_fit; RooFitResult *lau5_fit;
  RooFitResult *bern1_fit; RooFitResult *bern2_fit; RooFitResult *bern3_fit; RooFitResult *bern4_fit; RooFitResult *bern5_fit;
  CMS_hzg_mass.setBins(280);
  // RooFitResult *ggau_fit;
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
    pow1_fit = gauxpow1.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  pow3_fit = gauxpow3.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  pow5_fit = gauxpow5.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  exp1_fit = gauxexp1.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  exp3_fit = gauxexp3.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  exp5_fit = gauxexp5.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  lau1_fit = gauxlau1.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  lau2_fit = gauxlau2.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  lau3_fit = gauxlau3.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  lau4_fit = gauxlau4.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  lau5_fit = gauxlau5.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  // // // // // // bern1_fit = bern1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  bern2_fit = bern2.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(1)); //FIXME
	  bern3_fit = bern3.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  bern4_fit = bern4.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
	  bern5_fit = bern5.fitTo(datahist,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
    // ggau_fit = ggau.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));

  storedPdfs.add(gauxexp1);
  storedPdfs.add(gauxexp3);
  storedPdfs.add(gauxexp5);
  storedPdfs.add(gauxpow1);
  storedPdfs.add(gauxpow3);
  storedPdfs.add(gauxpow5);
  storedPdfs.add(gauxlau1);
  storedPdfs.add(gauxlau2);
  storedPdfs.add(gauxlau3);
  storedPdfs.add(gauxlau4);
  storedPdfs.add(gauxlau5);
  storedPdfs.add(bern2);
  storedPdfs.add(bern3);
  storedPdfs.add(bern4);
  storedPdfs.add(bern5);
  // storedPdfs.add(ggau);
  
  RooWorkspace *ws =  new RooWorkspace();
  RooCategory catIndex("catIndex","c");
  RooMultiPdf *pdf = new RooMultiPdf("CMS_hzg_bkgshape","all pdfs",catIndex,storedPdfs);
  RooRealVar nBackground("CMS_hzg_bkgshape_norm","nbkg",data.sumEntries(),0,10E8);
  ws->SetName("multipdf");
  ws->import(*pdf);
  ws->import(nBackground);
  ws->import(datahist);
  ws->import(catIndex);
  ws->import(data);
  fout->cd();
  ws->Write();
  bern2_fit->Write(); bern3_fit->Write(); bern4_fit->Write(); bern5_fit->Write();
  exp1_fit->Write(); exp3_fit->Write(); exp5_fit->Write();
  pow1_fit->Write(); pow3_fit->Write(); pow5_fit->Write();
  lau1_fit->Write(); lau2_fit->Write(); lau3_fit->Write();lau4_fit->Write();lau5_fit->Write();
  // ggau_fit->Write();

  fout->Close();
  RooPlot* xframe1  = CMS_hzg_mass.frame() ;
  CMS_hzg_mass.setRange("blind1",100,120) ;
  CMS_hzg_mass.setRange("blind2",130,170);
  data.plotOn(xframe1,Binning(70),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(70),CutRange("blind2")) ;
  gauxpow1.plotOn(xframe1,RooFit::Name("gauxpow1"),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashed));
  gauxpow3.plotOn(xframe1,RooFit::Name("gauxpow3"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashed));
  gauxpow5.plotOn(xframe1,RooFit::Name("gauxpow5"),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashed));
  gauxexp1.plotOn(xframe1,RooFit::Name("gauxexp1"),LineColor(TColor::GetColor("#FF8589")),LineStyle(kDotted));
  gauxexp3.plotOn(xframe1,RooFit::Name("gauxexp3"),LineColor(TColor::GetColor("#81C4EF")),LineStyle(kDotted));
  gauxexp5.plotOn(xframe1,RooFit::Name("gauxexp5"),LineColor(TColor::GetColor("#FFD35C")),LineStyle(kDotted));
  gauxlau1.plotOn(xframe1,RooFit::Name("gauxlau1"),LineColor(TColor::GetColor("#A30005")),LineStyle(kDashDotted));
  gauxlau2.plotOn(xframe1,RooFit::Name("gauxlau2"),LineColor(TColor::GetColor("#136090")),LineStyle(kDashDotted));
  gauxlau3.plotOn(xframe1,RooFit::Name("gauxlau3"),LineColor(TColor::GetColor("#E0A500")),LineStyle(kDashDotted));
  gauxlau4.plotOn(xframe1,RooFit::Name("gauxlau4"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashDotted));
  gauxlau5.plotOn(xframe1,RooFit::Name("gauxlau5"),LineColor(TColor::GetColor("#6A4C93")),LineStyle(kDashDotted));
  // bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(TColor::GetColor("#FF595E")));
  bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(TColor::GetColor("#FFCA3A")));
  bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(TColor::GetColor("#8AC926")));
  bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(TColor::GetColor("#1982C4")));
  bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(TColor::GetColor("#6A4C93")));
  // ggau.plotOn(xframe1, RooFit::Name("ggau"),LineColor(kGray+3));

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
  leg4->AddEntry(xframe1->findObject("gauxlau4"), "lau4", "l");
  leg4->AddEntry(xframe1->findObject("gauxlau5"), "lau5", "l");
  // leg4->AddEntry(xframe1->findObject("bern1"), "bern1", "l");
  leg4->AddEntry(xframe1->findObject("bern2"), "bern2", "l");
  leg4->AddEntry(xframe1->findObject("bern3"), "bern3", "l");
  leg4->AddEntry(xframe1->findObject("bern4"), "bern4", "l");
  leg4->AddEntry(xframe1->findObject("bern5"), "bern5", "l");
  // leg4->AddEntry(xframe1->findObject("ggau"), "ggau", "l");
    leg4->AddEntry(xframe1->findObject("data"), "data", "lep");

  leg4->Draw("same");
  gPad->Print(Form("cat%d_turn_lau.pdf",cat));
  cout<<"=======>Category"<<cat<<endl;
  
 cout<<"Bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();
 cout<<"Bern3:"<<bern3_fit->status()<<endl;bern3_fit->Print();
 cout<<"Bern4:"<<bern4_fit->status()<<endl;bern4_fit->Print();
 cout<<"Bern5:"<<bern5_fit->status()<<endl;bern5_fit->Print();
 cout<<"Pow1:"<<pow1_fit->status()<<endl;pow1_fit->Print();
 cout<<"Pow3:"<<pow3_fit->status()<<endl;pow3_fit->Print();
 cout<<"Pow5:"<<pow5_fit->status()<<endl;pow5_fit->Print();
 cout<<"Lau1:"<<lau1_fit->status()<<endl;lau1_fit->Print();
 cout<<"Lau2:"<<lau2_fit->status()<<endl;lau2_fit->Print();
 cout<<"Lau3:"<<lau3_fit->status()<<endl;lau3_fit->Print();
 cout<<"Lau4:"<<lau4_fit->status()<<endl;lau4_fit->Print();
 cout<<"Lau5:"<<lau5_fit->status()<<endl;lau5_fit->Print();
 cout<<"Exp1:"<<exp1_fit->status()<<endl;exp1_fit->Print();
 cout<<"Exp3:"<<exp3_fit->status()<<endl;exp3_fit->Print();
 cout<<"Exp5:"<<exp5_fit->status()<<endl;exp5_fit->Print();
//  cout<<"ggau:"<<ggau_fit->status()<<endl;ggau_fit->Print();
  // if(bern2_fit->status()!=0){cout<<"Bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();}
//   if(bern3_fit->status()!=0){cout<<"Bern3:"<<bern3_fit->status()<<endl;bern3_fit->Print();}
//   if(bern4_fit->status()!=0){cout<<"Bern4:"<<bern4_fit->status()<<endl;bern4_fit->Print();}
//   if(bern5_fit->status()!=0){cout<<"Bern5:"<<bern5_fit->status()<<endl;bern5_fit->Print();}
  // if(pow1_fit->status()!=0){cout<<"Pow1:"<<pow1_fit->status()<<endl;pow1_fit->Print();}
//   if(pow3_fit->status()!=0){cout<<"Pow3:"<<pow3_fit->status()<<endl;pow3_fit->Print();}
//   if(pow5_fit->status()!=0){cout<<"Pow5:"<<pow5_fit->status()<<endl;pow5_fit->Print();}
//  if(lau1_fit->status()!=0){cout<<"Lau1:"<<lau1_fit->status()<<endl;lau1_fit->Print();}
//   if(lau2_fit->status()!=0){cout<<"Lau2:"<<lau2_fit->status()<<endl;lau2_fit->Print();}
//   if(lau3_fit->status()!=0){cout<<"Lau3:"<<lau3_fit->status()<<endl;lau3_fit->Print();}
//   if(exp1_fit->status()!=0){cout<<"Exp1:"<<exp1_fit->status()<<endl;exp1_fit->Print();}
//   if(exp3_fit->status()!=0){cout<<"Exp3:"<<exp3_fit->status()<<endl;exp3_fit->Print();}
//   if(exp5_fit->status()!=0){cout<<"Exp5:"<<exp5_fit->status()<<endl;exp5_fit->Print();}
}
