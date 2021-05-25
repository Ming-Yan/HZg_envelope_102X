
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
void turnon_chi2(int cat,int nbind,int nbinu)
{
  // gSystem->Load("RooStepBernstein_cxx.so");
  // gSystem->Load("RooGaussStepBernstein_cxx.so");
  //gSystem->Load("/afs/cern.ch/work/m/milee/MYcode/limit/ModGaus11_cxx.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  // gSystem->Load("ModGaus11_cxx.so");
  TTree* DataTree1 = makeTTree(Form("bkg/ptwei_turnon/HiggsMass_ptwei_bkg_ele_mu_cat%d_2020_full.txt",cat));
  int totalev = DataTree1->GetEntriesFast();
  RooRealVar CMS_hzg_mass("CMS_hzg_mass", "CMS_hzg_mass", nbind, nbinu, "GeV") ;
  CMS_hzg_mass.setRange("window",nbind, nbinu);
  CMS_hzg_mass.setBins(4*(nbinu-nbind));
  RooDataSet data(Form("data_obs_ele_mu_cat%d_2020_13TeV",cat), " ", RooArgSet(CMS_hzg_mass), Import(*DataTree1));
  RooDataHist datahist(Form("CMS_hzg_datahist_ele_mu_cat%d_2020_13TeV",cat),Form("CMS_hzg_datahist_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,data);
  RooArgList storedPdfs(Form("CMS_hzg_datapdf_ele_mu_cat%d_2020_13TeV",cat));

  TFile *fout  = new TFile(Form("m%d_%d_turnon_cat%d_fix.root",nbind,nbinu,cat),"recreate");


  //generetic Bernstein polynomials
  RooRealVar mean(Form("mean_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("mean_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0);
  RooRealVar sigma_b1(Form("sigma_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  RooRealVar sigma_b2(Form("sigma_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  RooRealVar sigma_b3(Form("sigma_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  RooRealVar sigma_b4(Form("sigma_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  // RooRealVar sigma_b5(Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,8.);
  //cat4
  RooRealVar sigma_b5(Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.78,2.,8.);
  RooRealVar step_b1(Form("step_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,100,110);
  RooRealVar step_b2(Form("step_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,100,110);
  RooRealVar step_b3(Form("step_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,100,110);
  RooRealVar step_b4(Form("step_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105,100,110);
  RooRealVar step_b5(Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105,100,110);
 //cat4
  // RooRealVar step_b5(Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.58,100,110);
  RooRealVar p0(Form("p0_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p0_env_pdf_ele_mu_cat%d_2020_13TeV",cat),15);
  RooRealVar b1p1(Form("b1p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b1p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25);
  RooRealVar b2p1(Form("b2p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b2p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25);
  RooRealVar b2p2(Form("b2p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b2p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.,-25.,25);
  RooRealVar b3p1(Form("b3p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25);
  RooRealVar b3p2(Form("b3p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25);
  RooRealVar b3p3(Form("b3p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25.);
  RooRealVar b4p1(Form("b4p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25.);
  RooRealVar b4p2(Form("b4p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,-25.,25.);
  RooRealVar b4p3(Form("b4p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.1,-25.,25.);//untag
  RooRealVar b4p4(Form("b4p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.1,-25.,25);//VBF&lepton
  RooRealVar b5p1(Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.5,-15.,15.);//untag
  RooRealVar b5p2(Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.67,-15.,15.);
  RooRealVar b5p3(Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.24,-15.,15.);
  RooRealVar b5p4(Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.18,-15.,15.);
  RooRealVar b5p5(Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.17,-15.,15.);
  //cat4
  //  RooRealVar b5p1(Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.5663,-5.,5.);//untag
  // RooRealVar b5p2(Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.536,-5.,5.);
  // RooRealVar b5p3(Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.3253e-01,-15.,15.);
  // RooRealVar b5p4(Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.3505e-01,-15.,15.);
  // RooRealVar b5p5(Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.9944e-01,-15.,15.);
  RooGaussStepBernstein bern1(Form("bern1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_b1,step_b1, RooArgList(p0,b1p1));
  RooGaussStepBernstein bern2(Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_b2,step_b2, RooArgList(p0,b2p1,b2p2));
  RooGaussStepBernstein bern3(Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_b3,step_b3, RooArgList(p0,b3p1,b3p2,b3p3));
  RooGaussStepBernstein bern4(Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_b4,step_b4, RooArgList(p0,b4p1,b4p2,b4p3,b4p4));
  RooGaussStepBernstein bern5(Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_b5,step_b5, RooArgList(p0,b5p1,b5p2,b5p3,b5p4,b5p5));
  //testing with generic power law   
  //cat1
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.7,5.,8.);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.32,106,109);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.0,-8,-2.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.6576e-04,0.0,0.5);
  //cat2
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.0783,3.,10);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.5,107,110);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.1652,-15,-5.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.002 ,0.0,1.);
  //cat3
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4,3.,8.);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.6,106,107);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.5,-15,-5.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.02 ,0.0,1.);
  //cat4
  //  RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.);
  //  RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.,106,108);
  //  RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3453,-7,-5.);
  //  RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.7131e-05 ,0.0,1.);
  //cat501
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.,0.1,10);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,106,109);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.6,-15,-5.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0002 ,0.0,1.);
  // //cat502
  RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.,0.1,10);
  RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.2,104,110);
  RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.6,-8.,-5.);
  RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0006,0,1.);
  //cat503-100-170
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,0.1,5.);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.2,-8,-1.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.001,0.000001,1.);
  //105-170
  /*RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.8,1.,6.);
  RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104,105,1070);
  RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.9143,-8,-1.);
  RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.4993e-10,0.,1.);*/
  //cat6789
  // RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.95,1.,8.);
  // RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.7,104,106);
  // RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.,-8,-1.);
  // RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.001,0.000001,1.);
  //cat1
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.7641,5.,8.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.1,107,109);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4478,-8,-5.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.6823e-05,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.2,-7.,-2.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 2.8854e-11,0,1e-7);
  //cat2
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.6482,4.,7.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.91,107,110);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.4463,-9,-5);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.8038e-02,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7,-9.,-5.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.6736e-04,0,0.01);
  //cat3
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.6,106,107);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), -7.0864,-9,-6.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.8364,-6.,-3.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.1789e-14,0,0.001);
  //cat3-105-170
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.5483,3.,6.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.6,106,107);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), -7.4268,-9,-6.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9997e-01,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.8194,-7.,-3.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.2848e-16,0,0.000001);
  //cat4
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.5,3.,8.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,105,108);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3671,-10,-5.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.6543e-01,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.267,-8.,-2.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.1001e-06,0,1.);
  //105-170
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.0707,3.,5.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.87,106.,107);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4501,-8,-5.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.86,0.7,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.1,-10.,-3.5);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 6.2180e-02,0,1e-2);
  //cat501
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5779,7.,10.);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.45,108,109);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.7407,-10,-5.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.9305e-02,0.,0.5);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.6765,-8.,-4.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.4955e-08,0,0.0001);
  //cat502
  RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.,3.,10);
  RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105,100,110);
  RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.5,-10.,-5.);
  RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.7802e-05,0.,1.);
  RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.9,-6.,-1.);
  RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.2618e-01,0,1.);
  //cat503-100-170
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5,-8,-3.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4,-5.,-1.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.01,0,1.);
  //102-170
  // RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,7);
  // RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,106);
  // RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5,-8,-3.);
  // RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.3,0.,1.);
  // RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.9,-5.,-1.);
  // RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.01,0,1.);
  //cat6789
//  RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.5,3.,10.);
//   RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.39,104,107);
//   RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.2,-10.,-1.);
//   RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9959e-01,0.,1.);
//   RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.0,-10.,-1.);
//   RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0,.1);
  //cat1
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.6,5.,7.);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.42,107,109);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.2664,-7,-4);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 2.6960e-06,0.,.01);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.7909,-12.,-8.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.3933e-04,0,.1);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.0477,-8.,-5.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.5285e-02,0.,0.1);
  //cat2
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.3026,3.,8.);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.64,107,109);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.3134,-10,-5);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.e-01,0.,1.);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.0937,-10.,-2.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.0955e-01,0,1.);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.2482,-10.,-1.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.1590e-03,0.001,1.);
  //cat3
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.5,3.,6.5);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.6,106,107.5);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-8.4010,-11,-5);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.4091e-01,0.,1.);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.6841,-8.,-3.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2316e-04,0.,1.);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.9759,-8.,-3.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.0571e-06,0.0,1.);
  //cat4
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.0999,3.,5.);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.21,106,108);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.8,-11,-8);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.9425e-06,0.,0.001);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4921,-8.,-5.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.751e-02 ,0.,0.1);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3125,-8.,-5.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8914e-03,0.0,0.05);
  //cat501
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.4228,7.,12);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.75,107,109);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.6032,-8,-5);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.6078e-01,0.,1.);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4915,-8.,-2.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.4993e-03,0,1.);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4064,-8.,-2.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.1557e-02 ,0.001,1.);
  //cat502
  RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.4188,6.,8.5);
  RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.34,105,108);
  RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.1543,-8,-5);
  RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.6163e-01,0.0,1.);
  RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.2342,-8.,-4.);
  RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.1945e-01,0.1,1.);
  RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.,-8.,-4.);
  RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.7946e-02,0.001,.5);
  //cat503
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,10);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6,-12,-5);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.07,0.,1.);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.,-10.,-2.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.03,0,1.);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4,-10.,-1.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0001,0.00001,1.);
  // RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,7);
  // RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,107);
  // RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4019,-12,-5);
  // RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.1582e-02,0.,1.);
  // RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-8.6290,-10.,-2.);
  // RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.5799e-04 ,0,0.1);
  // RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.8322,-10.,-1.);
  // RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.2964e-04,0.00,.1);
  //cat6789
//  RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.95,3.,10);
//   RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.9,104,107);
//   RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.2,-12,-4);
//   RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0.,1.);
//   RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.,-10.,0.);
//   RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.37,0,1.);
//   RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.8,-10.,0.);
//   RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0.0,1.);

  RooGenericPdf step_pow1(Form("step_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(CMS_hzg_mass,turnon_pow1,p1_pow1,cp1_pow1));//step*(ax^b)
  RooGenericPdf step_pow3(Form("step_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(CMS_hzg_mass,turnon_pow3,p1_pow3,cp1_pow3,p3_pow3,cp3_pow3));//step*(ax^b+cx^d)
  RooGenericPdf step_pow5(Form("step_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(CMS_hzg_mass,turnon_pow5,p1_pow5,cp1_pow5,p3_pow5,cp3_pow5,p5_pow5,cp5_pow5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_pow1(Form("gau_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_pow1);
  RooGaussModel gau_pow3(Form("gau_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_pow3);
  RooGaussModel gau_pow5(Form("gau_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_pow5);
  RooFFTConvPdf gauxpow1(Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_pow1,gau_pow1);
  RooFFTConvPdf gauxpow3(Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_pow3,gau_pow3);
  RooFFTConvPdf gauxpow5(Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_pow5,gau_pow5);
  gauxpow1.setBufferFraction(0.25);
  gauxpow3.setBufferFraction(0.25);
  gauxpow5.setBufferFraction(0.25);
  //testing with generic Laurent
  //cat1
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.8148,3.,10);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,105,109);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.7224e-14,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.1208e-01,0,1.);
  //cat2
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.0908,4.,6.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.1,106,109);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.5039e-06,0.0,0.1);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.3050e-01,0,1.);
  //cat3
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,10.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,105.,109);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.6051e-13 ,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.0265e-01,0,1.);
  //cat4
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.6496,2.,6.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.05,105,109);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.e-12,0.0,1e-10);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.7 ,0.6,1.);
  //cat501
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.000000001,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99999,0,1.);
  //cat502
  RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,100,110);
  RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.000000001,0.0,1.);//untag
  RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99999,0,1.);
  //cat503
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.9,3.,5.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.8,104,106);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.3937e-01,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.9829e-06,0,1.);
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.9,3.,7.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.5,104,107);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5474e-01,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.6236e-08,0,.1);
  //cat6789
  // RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.95,3.,6.);
  // RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.9,103,107);
  // RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.7933e-01,0.0,1.);//untag
  // RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.3706e-03,0,1.);
  //cat1
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000000000000000000000000000000000001,0.,1.);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.9999999999999999999999999999999,0.,1.05);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,1.);
  //cat2
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.0,4.,5.5);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2503e-06,0.,1.);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.7549e-01,0.1,1.);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.1416e-09,0,1.);
  //cat3
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.,106.,108.);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.0493e-07,0.,0.1);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.8448e-01 ,0.5,1.0);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.4361e-08,0,0.1);
  //cat4
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.1145,2.,6);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.89,105,109);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2209e-05,0.,0.01);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.2780e-01,0.5,1.0);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2371e-07,0,0.01);
  //cat501
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7,3.,10);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.1098e-03,0.01,1.);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0.,1.);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,0.00001);
  //cat502
  RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,100,110);
  RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.00000000001,0.,1.);
  RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0.,1.);
  RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,1.);
  // //cat503
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,6.);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,106);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.7839e-06,0.,1.);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.95,0.,1.05);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.1,0,1.);
  //cat6789
  // RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.95,3.,6.);
  // RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.9,104,107);
  // RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.7839e-06,0.,1.);
  // RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.95,0.,1.05);
  // RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.1,0,1.);
  //cat1
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000000001,0.,1.);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0.,1.);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,1.);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  //cat2
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.8,3.,6.);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.78,107,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.e-09,0.,.001);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.e-08,0.,.1);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7e-11,0,.0001);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.1e-01,0.,1.);
  // cat3
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,1.,6.);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,104,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8872e-15,0.,.001);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.2066e-13,0.,.1);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5240e-17,0,.0001);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9935e-01,0.9,1.);
  //105-170
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.4211,1.,6.);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.39,104,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.4094e-12,0.,1e-6);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.6946e-07,0.,.01);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.7928e-13,0,1e-9);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.0139e-01,0.5,1.);
  //cat4
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.53,2.,6.);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.2,105,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.3514e-10,0.,.001);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5466e-08,0.,.1);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.7939e-11,0,.0001);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.0697e-01,0.9,1.);
  //cat501
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.7,3.,10);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,105,109);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.2769e-07,0.,0.00001);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.4839e-01,0.,1.);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.4089e-07,0,0.00001);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  // cat502
  RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.,3.,10);
  RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.34,100,110);
  RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.7439e-12,0.,1e-8);
  RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5179e-03,0.,0.01);
  RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.4468e-08,0,0.000001);
  RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99,0,1.);
  //cat503
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,8);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0.,1.);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.06,0.,1.);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.000001,0,1.);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.02,0,1.);
  //cat6789
  // RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,6.);
  // RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.,104,107);
  // RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.0862e-02,0.,1.);
  // RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.1937e-02,0.,1.);
  // RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.7030e-03,0,1.);
  // RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.5681e-01,0,1.);
  
  //cat1
  //  RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.3881,3.,10);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9874e-08,0.,0.0000001);  
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9867e-07,0.,0.0001);   
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.9584e-09,0,0.00000001);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.97335,0.5,1.); 
  // RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.3252e-11 ,0,0.000000001); 
  //cat2
//   RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.4050,3.,6.);
//   RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.75,107,109);
//   RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.0776e-08,0.,.0001);
//   RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 3.2590e-06,0.,.01);
//   RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.7575e-12,0,.0001);
//   RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.0505e-01,0.,1.);
//  RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.7040e-13 ,0,.00001);

  // cat3
  // RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.4105,1.,6.);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.49,104,109);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.6620e-09,0.,1e-5);
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.3996e-08,0.,1e-3);
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.0450e-11,0,1e-7);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.4883e-01,0.,1.);
  // RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8935e-14,0,1e-10);
  //cat4
  // RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.53,2.,6.);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.2,105,109);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.3514e-05,0.,.001);
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5466e-02,0.,.1);
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.7939e-17,0,.0001);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99,0.9,1.);
  // RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,1.);
  //cat501
  // RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000000001,0.,.0001);
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0.,0.01);
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,0.0000001);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  // RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,.000001);
  // cat502
  RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,5.5);
  RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.1,105,108);
 RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.4302e-07,0.,.001);
  RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5957e-03,0.,.1);
  RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.9832e-10,0,.000001);
  RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9291e-01,0.99,1.);
  RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,.000000000001);
  //cat503
  // RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,8);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0.,1.);
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.06,0.,1.);
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.000001,0,1.);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.02,0,1.);
  //   RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,1.);
    //cat6789
  // RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,6.);
  // RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.,104,107);
  // RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.0862e-02,0.,1.);
  // RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.1937e-02,0.,1.);
  // RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.7030e-03,0,1.);
  // RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.5681e-01,0,1.);
  // RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,1.);
  //cat1
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.6367,5.,7);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.6749e-11,0.,0.0000001);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.1346e-05,0.,0.0001);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.3162e-14,0,0.00000001);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.6,0,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.6445e-12,0,1e-05);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999,0,1.);
    //cat2
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.8,3.,6.);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.78,107,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.e-09,0.,.001);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.e-08,0.,.1);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7e-11,0,.0001);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.1e-01,0.,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.00000,0,0.00001);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.00000,0,0.00001);
  // cat3
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,6.);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,104,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.e-09,0.,.001);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.e-08,0.,.1);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7e-11,0,.0001);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.1e-01,0.,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0,1.);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0,1.);
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.5663,2.,6.);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.86,104,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.0048e-06,0.,1e-4);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.8894e-02,0.,.1);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.2703e-08,0,1e-5);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.8445e-01,0.,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.9377e-11,0,1e-6);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.5449e-01,0,1.);
  //cat4
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.3253,2.,6.);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.41,105,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.8677e-05,0.,.1);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.9106e-01,0.,.9);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.7235e-09,0,.0001);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.1470e-01 ,0.1,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.5457e-09,0,1e-4);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.5973e-02,0.0,0.5);
  //cat501
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6,3.,10);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108,107,109);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000000001,0.,0.001);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0000001,0.,0.01);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.0,0,0.01);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.4,0,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,0.0000001);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  // cat502
  RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5,3.,10);
  RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.1,105,108);
  RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5897e-06,0.,0.001);
  RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9863e-03,0.,0.1);
  RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9971e-07,0,0.0001);
  RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.0000e-10,0,0.00000001);
  RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.01,0,.1);
  RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  //cat503
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,8);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0.,1.);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.06,0.,1.);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.000001,0,1.);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.02,0,1.);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,1.);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.999999999,0,1.);
  //cat6789
  // RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.,3.,6.);
  // RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.,104,107);
  // RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.1069e-01 ,0.,1.);
  // RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.1937e-02,0.,1.);
  // RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.3526e-01,0,1.);
  // RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.8827e-04,0,.1);
  // RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.4479e-08,0,.1);
  // RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.8845e-02,0,1.);

  RooGenericPdf step_lau1(Form("step_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(CMS_hzg_mass,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
  RooGenericPdf step_lau2(Form("step_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(CMS_hzg_mass,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
  RooGenericPdf step_lau3(Form("step_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(CMS_hzg_mass,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
  RooGenericPdf step_lau4(Form("step_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2))", RooArgList(CMS_hzg_mass,turnon_lau4,cl1_lau4,cl2_lau4,cl3_lau4,cl4_lau4,cl5_lau4));//step*(ax^b+cx^d)
  RooGenericPdf step_lau5(Form("step_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2)+@7*(@0)^(-7))", RooArgList(CMS_hzg_mass,turnon_lau5,cl1_lau5,cl2_lau5,cl3_lau5,cl4_lau5,cl5_lau5,cl6_lau5));//step*(ax^b+cx^d)

  RooGaussModel gau_lau1(Form("gau_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_lau1);
  RooGaussModel gau_lau2(Form("gau_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_lau2);
  RooGaussModel gau_lau3(Form("gau_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_lau3);
  RooGaussModel gau_lau4(Form("gau_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_lau4);
  RooGaussModel gau_lau5(Form("gau_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_lau5);

  RooFFTConvPdf gauxlau1(Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_lau1,gau_lau1);
  RooFFTConvPdf gauxlau2(Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_lau2,gau_lau2);
  RooFFTConvPdf gauxlau3(Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_lau3,gau_lau3);
  RooFFTConvPdf gauxlau4(Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_lau4,gau_lau4);
  RooFFTConvPdf gauxlau5(Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_lau5,gau_lau5);
  gauxlau1.setBufferFraction(0.25);
  gauxlau2.setBufferFraction(0.25);
  gauxlau3.setBufferFraction(0.25);
  gauxlau4.setBufferFraction(0.25);
  gauxlau5.setBufferFraction(0.25);

  // testing with generic exponential
  //cat1
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.3460  ,5.,8);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 107.5,107,109);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.3265e-02,-0.7,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0,1.);
  //cat2
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.9697,4.,6.);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.07,107.,110);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.2395e-02,-0.2,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.1021e-01,0,1.);
  //cat3
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,106.,108.);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.05,-0.5,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0,0.1);
  //cat4
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.56,3.,4.);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.97,106.,108);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4e-02,-0.1,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.01,0,1.);
  //cat501
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5,3.,10);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 108,107,109);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.03,-0.7,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0,1.);
  //cat502
  RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5,3.,10);
  RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,100,110);
  RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.04,-0.7,0.);
  RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0,1.);
  //cat503
  //RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4,3.,10);
  //RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 104.6,104,105);
  //RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.03,-0.7,0.);
  //RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0,1.);
  //cat6789
  // RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.0,5.,6.);
  // RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 106.3,105,106.5);
  // RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.0228,-0.7,0.);
  // RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.023,0,1.);

  //cat1
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.2293,6.,7);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.85,107,109);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.0820e-02,-0.7,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.3882e-02,0.,1.);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), -4.2712e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.2016e-02,0,1.);
  //cat2
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.02,3.,6.);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.41,100,110);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.6444e-02,-0.5,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.7480e-01,0.,1.);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.1990e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.1965e-02,0,.5);
  //cat3
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.5);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,105,108.);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.4503e-02,-0.5,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.6120e-01,0.,1.);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.8149e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8388e-07,0,.1);
  //m105-170
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.3635,3.,5.5);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.97,106,108.);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.6239e-02,-0.7,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8116e-03,0.,1.);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.2833e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9987e-04 ,0,1e-2);
  //cat4
  //  RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,3.,5.);
  //  RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.5,105,108);
  //  RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.4437e-02,-0.5,0.);
  //  RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2335e-02,0.,1.);
  //  RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.5070e-02,-0.5,0.);
  //  RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.9216e-01,0,1.);
  //cat501
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.,7.,12);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.5,107,109);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.2614e-02,-0.5,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.6306e-09,0.,0.0001);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.7779e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.8390e-05,0,0.01);
  //cat502
  RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.2184,5.,8);
  RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.75,104.6,110);
  RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), -4.2557e-02,-0.5,0.);
  RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.01,0.00001,1.);
  RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.4618e-02 ,-1.,0.);
  RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.9999e-01,0.01,1.);
  //cat503
  //RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.832,2.,6.);
  //RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
   //RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.025,-0.1,0.);
   //RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.17,0.,1.);
   //RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.024,-0.1,0.);
   //RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.1120e-06,0,1.);
  //cat6789
  // RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.4,4.,7.);
  // RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.2,104,107);
  // RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.5777e-02,-0.1,0.);
  // RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.6637e-01,0.,1.);
  // RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-2e-02,-0.5,0.);
  // RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.e-03,0,0.1);
  
  //cat1
//   RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.22,6,7);
//   RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.8,107,108);
//  RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.0820e-02,-1.,0.);
//   RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), 4.6641e-02,0.,1.);
//   RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.2712e-02,-0.5,0.);
//   RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.7630e-02,0,0.1);
//   RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.5734e-03,-0.1,0.);
//   RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.0000e-05,0,.001);
  //cat2
  // RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.0658,4.,6.);
  // RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),108.3,107,109);
  // RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.4155e-02,-0.5,0.);
  // RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.7451e-01,0.,1.);
  // RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.5821e-02,-0.5,0.);
  // RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.3393e-07,0,.001);
  // RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-2.4270e-02,-0.5,0.);
  // RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.8593e-04,0,.1);
  //cat3
  // RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.85,3.,6);
  // RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,105,118);
  // RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.054,-0.5,0.);
  // RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.5202e-01,0.,1.);
  // RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.055,-0.5,0.);
  // RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.9097e-01,0,1.);
  // RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-0.059,-0.5,0.);
  // RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.9429e-09,0,1.);
  //cat4
  //  RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.12,3.,6);
  //  RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.21,106.5,108);
  //  RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.1816e-02,-0.5,0.);
  //  RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0.0,1.);
  //  RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.3062e-02,-0.5,0.);
  //  RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0,.1);
  //  RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.2008e-04,-0.01,0.);
  //  RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.5533e-05,0,.001);
   
  //cat501
  // RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.9229,7.,10);
  // RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.5,107,109);
  // RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.7776e-02,-0.5,0.);
  // RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.9302e-07,0.,0.0001);
  // RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.7676e-02,-0.5,0.);
  // RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8750e-04,0,0.01);
  // RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.7448e-04 ,-0.001,0.);
  // RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.,0,1e-8);
  //cat502
  RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.8258,5.,7.5);
  RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.31,104,106.5);
  RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.3962e-02 ,-0.5,0.);
  RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.2914e-02,0.,0.5);
  RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.4028e-02 ,-0.5,0.);
  RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.2306e-02,0,0.5);
  RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4e-02,-0.5,0.);
  RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_ex  p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),7.6517e-09,0,1e-4);
  //cat503
   /*RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.8366,3.,5.);
  RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.6,104,105);
  RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-1.5800e-01,-0.5,0.);
  RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.9251e-04,0.,1.);
  RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-2.4977e-02,-0.5,0.);
  RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.6553e-03,0,1.);
  RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.7189e-03,-0.5,0.);
  RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.6356e-07 ,0,1.);*/
  //cat6789
  // RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.56,3.,6.);
  // RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.3,104,107);
  // RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-2.7150e-02,-0.6,0.);
  // RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.4506e-01,0.,1.);
  // RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-2.4988e-02,-0.5,0.);
  // RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.7690e-06,0,.1);
  // RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-7.9924e-06,-0.5,0.);
  // RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.1386e-03,0,.1);
  RooGenericPdf step_exp1(Form("step_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(CMS_hzg_mass,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
  RooGenericPdf step_exp3(Form("step_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(CMS_hzg_mass,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
  RooGenericPdf step_exp5(Form("step_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(CMS_hzg_mass,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_exp1(Form("gau_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_exp1);
  RooGaussModel gau_exp3(Form("gau_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_exp3);
  RooGaussModel gau_exp5(Form("gau_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,mean,sigma_exp5);
  RooFFTConvPdf gauxexp1(Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_exp1,gau_exp1);
  RooFFTConvPdf gauxexp3(Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_exp3,gau_exp3);
  RooFFTConvPdf gauxexp5(Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),CMS_hzg_mass,step_exp5,gau_exp5);
  gauxexp1.setBufferFraction(0.25);
  gauxexp3.setBufferFraction(0.25);
  gauxexp5.setBufferFraction(0.25);
  //Generic Gaussian
  // RooRealVar turnon_ggau(Form("turnon_ggau",Form("turnon_ggau",112,108,120);
  // RooRealVar turnon_gv0(Form("turnon_gv0",Form("turnon_gv0",1.6,0.,5.);
  // RooRealVar turnon_gv1(Form("turnon_gv1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_gv1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1e-10,0.,0.5);
  // RooRealVar turnon_gs0(Form("turnon_gs0",Form("turnon_gs0",1.5,1.,10);
  // RooRealVar turnon_gs1(Form("turnon_gs1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_gs1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),28,5.,50.);
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
  // // // // // // // bern1_fit = bern1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
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
  RooCategory catIndex(Form("CMS_hzg_pdfindex_ele_mu_cat%d_2020_13TeV",cat),Form("CMS_hzg_pdfindex_ele_mu_cat%d_2020_13TeV",cat));
  RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hzg_bkgshape_ele_mu_cat%d_2020_13TeV",cat),"all pdfs",catIndex,storedPdfs);
  RooRealVar nBackground(Form("CMS_hzg_bkgshape_ele_mu_cat%d_2020_13TeV_norm",cat),"nbkg",data.sumEntries(),0,10E8);
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
  xframe1->SetTitle("");
  xframe1->SetXTitle("m_{l^{+}l^{-}#gamma}");
  CMS_hzg_mass.setRange("blind1",nbind,120) ;
  CMS_hzg_mass.setRange("blind2",130,nbinu);
  data.plotOn(xframe1,Binning(nbinu-nbind),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(nbinu-nbind),CutRange("blind2")) ;
  gauxpow1.plotOn(xframe1,RooFit::Name(Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashed));
  gauxpow3.plotOn(xframe1,RooFit::Name(Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashed));
  gauxpow5.plotOn(xframe1,RooFit::Name(Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashed));
  gauxexp1.plotOn(xframe1,RooFit::Name(Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDotted));
  gauxexp3.plotOn(xframe1,RooFit::Name(Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#81C4EF")),LineStyle(kDotted));
  gauxexp5.plotOn(xframe1,RooFit::Name(Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FFD35C")),LineStyle(kDotted));
  gauxlau1.plotOn(xframe1,RooFit::Name(Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashDotted));
  gauxlau2.plotOn(xframe1,RooFit::Name(Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashDotted));
  gauxlau3.plotOn(xframe1,RooFit::Name(Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#8AC926")),LineStyle(kDashDotted));
  gauxlau4.plotOn(xframe1,RooFit::Name(Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashDotted));
  gauxlau5.plotOn(xframe1,RooFit::Name(Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#6A4C93")),LineStyle(kDashDotted));
  bern2.plotOn(xframe1, RooFit::Name(Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#FFCA3A")));
  bern3.plotOn(xframe1, RooFit::Name(Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#8AC926")));
  bern4.plotOn(xframe1, RooFit::Name(Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#1982C4")));
  bern5.plotOn(xframe1, RooFit::Name(Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)),LineColor(TColor::GetColor("#6A4C93")));
  // ggau.plotOn(xframe1, RooFit::Name("ggau"),LineColor(kGray+3));
  xframe1->SetMinimum(0.0001);
  xframe1->SetTitleSize(0.04);
  xframe1->SetLabelSize(0.035);
  xframe1->Draw();
  TLegend* leg4 = new TLegend(0.7,0.5,0.9,0.9);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  // leg4->SetNColumns(3);
  leg4->AddEntry(xframe1->findObject(Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxpow1", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxpow3", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxpow5", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxexp1", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxexp3", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxexp5", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxlau1", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxlau2", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxlau3", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxlau4", "l");
  leg4->AddEntry(xframe1->findObject(Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxlau5", "l");
  leg4->AddEntry(xframe1->findObject(Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxbern2", "l");
  leg4->AddEntry(xframe1->findObject(Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxbern3", "l");
  leg4->AddEntry(xframe1->findObject(Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxbern4", "l");
  leg4->AddEntry(xframe1->findObject(Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat)), "gauxbern5", "l");
  leg4->AddEntry(xframe1->findObject("ggau"), "ggau", "l");
  leg4->AddEntry(xframe1->findObject(Form("data",cat)), "data", "lep");

  leg4->Draw("same");
  gPad->Print(Form("m%d_%d_cat%d_turn_lau.pdf",nbind,nbinu,cat));
  cout<<"=======>Category"<<cat<<endl;
 
  cout<<"bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();
  cout<<"bern3:"<<bern3_fit->status()<<endl;bern3_fit->Print();
  cout<<"bern4:"<<bern4_fit->status()<<endl;bern4_fit->Print();
  cout<<"bern5:"<<bern5_fit->status()<<endl;bern5_fit->Print();
  cout<<"Pow1:"<<pow1_fit->status()<<endl;pow1_fit->Print();
  cout<<"Pow3:"<<pow3_fit->status()<<endl;pow3_fit->Print();
  cout<<"Pow5:"<<pow5_fit->status()<<endl;pow5_fit->Print();
  cout<<"Lau1:"<<lau1_fit->status()<<endl;lau1_fit->Print();
  cout<<"Lau2:"<<lau2_fit->status()<<endl;lau2_fit->Print();
  cout<<"Lau3:"<<lau3_fit->status()<<endl;lau3_fit->Print();
  cout<<"Lau4:"<<lau4_fit->status()<<" "<<lau4_fit->minNll()<<endl;lau4_fit->Print();
  cout<<"Lau5:"<<lau5_fit->status()<<" "<<lau5_fit->minNll  ()<<endl;lau5_fit->Print();
  cout<<"Exp1:"<<exp1_fit->status()<<endl;exp1_fit->Print();
  cout<<"Exp3:"<<exp3_fit->status()<<endl;exp3_fit->Print();
  cout<<"Exp5:"<<exp5_fit->status()<<endl;exp5_fit->Print();
  //  cout<<"ggau:"<<ggau_fit->status()<<endl;ggau_fit->Print();
  
}
