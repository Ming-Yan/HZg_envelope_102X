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

  TTree* DataTree1 = makeTTree(Form("HiggsMass_ptwei_bkg_ele_mu_cat%d_2020_full.txt",cat));
  int totalev = DataTree1->GetEntriesFast();
  RooRealVar mH("mH", "mH", 100, 180, "GeV") ;
  mH.setRange("window",100, 180);
  RooDataSet data(Form("data_obs_ele_mu_cat%d_2020",cat), " ", RooArgSet(mH), Import(*DataTree1));
  RooArgList storedPdfs("store");
  TFile *fout  = new TFile(Form("turnon_cat%d.root",cat),"recreate");
  //generetic Bernstein polynomials
  
  RooRealVar mean("mean","mean",0);
  RooRealVar sigma_b1("sigma_b1","sigma_b1",3,0.1,10);
  RooRealVar sigma_b2("sigma_b2","sigma_b2",3,0.1,10);
  RooRealVar sigma_b3("sigma_b3","sigma_b3",3,0.1,10);
  RooRealVar sigma_b4("sigma_b4","sigma_b4",3,0.1,10);
  RooRealVar sigma_b5("sigma_b5","sigma_b5",3,0.1,10);
  RooRealVar step_b1("step_b1","step_b1",105,90,120);
  RooRealVar step_b2("step_b2","step_b2",105,90,120);
  RooRealVar step_b3("step_b3","step_b3",105,90,120);
  RooRealVar step_b4("step_b4","step_b4",105,90,120);
  RooRealVar step_b5("step_b5","step_b5",105,90,120);
  RooRealVar p0("p0","p0",15);
  RooRealVar b1p1("b1p1","b1p1",0.3,-5.,5);
  RooRealVar b2p1("b2p1","b2p1",0.3,-5.,5);
  RooRealVar b3p1("b3p1","b3p1",0.3,-5.,5);
  RooRealVar b4p1("b4p1","b4p1",0.3,-5.,5);
  // RooRealVar b5p1("b5p1","b5p1",-0.1,-5.,5);//VBF&lepton
  RooRealVar b5p1("b5p1","b5p1",0.5,-5.,5);//untag
  RooRealVar b2p2("b2p2","b2p2",0.3,-5.,5);
  RooRealVar b3p2("b3p2","b3p2",0.3,-5.,5);
  RooRealVar b4p2("b4p2","b4p2",0.3,-5.,5);
  RooRealVar b5p2("b5p2","b5p2",0.1,-5.,5);
  RooRealVar b3p3("b3p3","b3p3",0.3,-5.,5);
  // RooRealVar b4p3("b4p3","b4p3",0.3,-5.,5);//VBF&lepton
  RooRealVar b4p3("b4p3","b4p3",-0.1,-5.,5);//untag
  RooRealVar b5p3("b5p3","b5p3",0.7,-5.,5);
  RooRealVar b3p4("b3p4","b3p4",0.3,-5.,5);
  // RooRealVar b4p4("b4p4","b4p4",0.3,-5.,5);//VBF&lepton
  RooRealVar b4p4("b4p4","b4p4",0.1,-5.,5);//VBF&lepton
  RooRealVar b5p4("b5p4","b5p4",0.5,-5.,5);
  RooRealVar b3p5("b3p5","b3p5",0.3,-5.,5);
  RooRealVar b4p5("b4p5","b4p5",0.3,-5.,5);
  RooRealVar b5p5("b5p5","b5p5",0.5,-5.,5);
  RooGaussStepBernstein bern1("bern1","bern1",mH,mean,sigma_b1,step_b1, RooArgList(p0,b1p1));
  RooGaussStepBernstein bern2("bern2","bern2",mH,mean,sigma_b2,step_b2, RooArgList(p0,b2p1,b2p2));
  RooGaussStepBernstein bern3("bern3","bern3",mH,mean,sigma_b3,step_b3, RooArgList(p0,b3p1,b3p2,b3p3));
  RooGaussStepBernstein bern4("bern4","bern4",mH,mean,sigma_b4,step_b4, RooArgList(p0,b4p1,b4p2,b4p3,b4p4));
  RooGaussStepBernstein bern5("bern5","bern5",mH,mean,sigma_b5,step_b5, RooArgList(p0,b5p1,b5p2,b5p3,b5p4,b5p5));
  //testing with generic power law 
  RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",5.,0.1,10);
  RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",107,90,120);
  RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",5,0.1,10);
  RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",107,100,120);
  RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",5,0.1,10);
  RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",107,90,120);
  RooRealVar p1_pow1("p1_pow1","p1_pow1",-6.5,-15,-5.);
  RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.3,0,1.);
  RooRealVar p1_pow3("p1_pow3","p1_pow3",-6.5,-15,-5.);
  RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",0.4,0.,1.);
  RooRealVar p3_pow3("p3_pow3","p3_pow3",-6,-10.,-2.);
  RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.01,0,1.);
  RooRealVar p1_pow5("p1_pow5","p1_pow5",-9,-15,-5);
  RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",0.7,0.,1.);
  RooRealVar p3_pow5("p3_pow5","p3_pow5",-8.,-10.,-2.);
  RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.1,0,1.);
  RooRealVar p5_pow5("p5_pow5","p5_pow5",-5.,-10.,-1.);
  RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.01,0.001,1.);
  RooGenericPdf step_pow1("step_pow1", "step_pow1", "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(mH,turnon_pow1,p1_pow1,cp1_pow1));//step*(ax^b)
  RooGenericPdf step_pow3("step_pow3", "step_pow3", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(mH,turnon_pow3,p1_pow3,cp1_pow3,p3_pow3,cp3_pow3));//step*(ax^b+cx^d)
  RooGenericPdf step_pow5("step_pow5", "step_pow5", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(mH,turnon_pow5,p1_pow5,cp1_pow5,p3_pow5,cp3_pow5,p5_pow5,cp5_pow5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_pow1("gau_pow1","gau_pow1",mH,mean,sigma_pow1);
  RooGaussModel gau_pow3("gau_pow3","gau_pow3",mH,mean,sigma_pow3);
  RooGaussModel gau_pow5("gau_pow5","gau_pow5",mH,mean,sigma_pow5);
  RooFFTConvPdf gauxpow1("gauxpow1","gauxpow1",mH,step_pow1,gau_pow1);
  RooFFTConvPdf gauxpow3("gauxpow3","gauxpow3",mH,step_pow3,gau_pow3);
  RooFFTConvPdf gauxpow5("gauxpow5","gauxpow5",mH,step_pow5,gau_pow5);
  
  //testing with generic Laurent
  RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",4,0.1,10);
  RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",105,100,120);
  RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",4,0.1,10);
  RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",105,100,120);
  RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",4,0.1,10);
  RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",105,100,120);
  // RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.25,0,1.);//VBF&lep
  RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.9,0.01,1.2);//untag
  RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.99,0,1.2);
  RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.93  ,0.,1.2);
  RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.99,0.,1.2);
  RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.01,0,1.2);
  RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.8,0.,1.2);
  RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.9,0.,1.2);
  RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.1,0,1.2);
  RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.99,0,1.2);
  RooGenericPdf step_lau1("step_lau1", "step_lau1", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(mH,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
  RooGenericPdf step_lau2("step_lau2", "step_lau2", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(mH,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
  RooGenericPdf step_lau3("step_lau3", "step_lau3", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(mH,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
  RooGaussModel gau_lau1("gau_lau1","gau_lau1",mH,mean,sigma_lau1);
  RooGaussModel gau_lau2("gau_lau2","gau_lau2",mH,mean,sigma_lau2);
  RooGaussModel gau_lau3("gau_lau3","gau_lau3",mH,mean,sigma_lau3);
  RooFFTConvPdf gauxlau1("gauxlau1","gauxlau1",mH,step_lau1,gau_lau1);
  RooFFTConvPdf gauxlau2("gauxlau2","gauxlau2",mH,step_lau2,gau_lau2);
  RooFFTConvPdf gauxlau3("gauxlau3","gauxlau3",mH,step_lau3,gau_lau3);


  // testing with generic exponential
  RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",5.5,0.1,10);
  RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",107,100,110);
  RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",5.5,0.1,10);
  RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",107,100,110);
  RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",5.5,0.1,10);
  RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",107,90,120);
  RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.05,-0.7,0.);
  RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.05,-0.5,0.);
  RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",0.4,0.,1.);
  RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.05,-0.5,0.);
  RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.5,0,1.);
  RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.15,-0.5,0.);
  RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",0.8,0.,1.);
  RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.09 ,-0.5,0.);
  RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.1,0,1.);
  RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.05,-0.5,0.);
  RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.1,0,1.);
  RooGenericPdf step_exp1("step_exp1", "step_exp1", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(mH,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
  RooGenericPdf step_exp3("step_exp3", "step_exp3", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(mH,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
  RooGenericPdf step_exp5("step_exp5", "step_exp5", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(mH,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_exp1("gau_exp1","gau_exp1",mH,mean,sigma_exp1);
  RooGaussModel gau_exp3("gau_exp3","gau_exp3",mH,mean,sigma_exp3);
  RooGaussModel gau_exp5("gau_exp5","gau_exp5",mH,mean,sigma_exp5);
  RooFFTConvPdf gauxexp1("gauxexp1","gauxexp1",mH,step_exp1,gau_exp1);
  RooFFTConvPdf gauxexp3("gauxexp3","gauxexp3",mH,step_exp3,gau_exp3);
  RooFFTConvPdf gauxexp5("gauxexp5","gauxexp5",mH,step_exp5,gau_exp5);
  
  
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
	// while (stat>1){
	//   if (ntries>=5) break;
	  pow1_fit = gauxpow1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = pow1_fit->status();
	//   minnll = pow1_fit->minNll();
	//   if (stat>1) params_pow1->assignValueOnly(pow1_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=5 ) break;
	  pow3_fit = gauxpow3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = pow3_fit->status();
	//   minnll = pow3_fit->minNll();
	//   if (stat>1) params_pow3->assignValueOnly(pow3_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=5) break;
	  pow5_fit = gauxpow5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = pow5_fit->status();
	//   minnll = pow5_fit->minNll();
	//   if (stat>1) params_pow5->assignValueOnly(pow5_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  exp1_fit = gauxexp1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = exp1_fit->status();
	//   minnll = exp1_fit->minNll();
	//   if (stat>1) params_exp1->assignValueOnly(exp1_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  exp3_fit = gauxexp3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = exp3_fit->status();
	//   minnll = exp3_fit->minNll();
	//   if (stat>1) params_exp3->assignValueOnly(exp3_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  exp5_fit = gauxexp5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = exp5_fit->status();
	//   minnll = exp5_fit->minNll();
	//   if (stat>1) params_exp5->assignValueOnly(exp5_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  lau1_fit = gauxlau1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = lau1_fit->status();
	//   minnll = lau1_fit->minNll();
	//   if (stat>1) params_lau1->assignValueOnly(lau1_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  lau2_fit = gauxlau2.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = lau2_fit->status();
	//   minnll = lau2_fit->minNll();
	//   if (stat>1) params_lau2->assignValueOnly(lau2_fit->randomizePars());
	//   ntries++; 
	// }
  // stat=1;//initialize 
  // ntries=0;//initialize
  // minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  lau3_fit = gauxlau3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = lau3_fit->status();
	//   minnll = lau3_fit->minNll();
	//   if (stat>1) params_lau3->assignValueOnly(lau3_fit->randomizePars());
	//   ntries++; 
	// }
  stat=1;//initialize 
  ntries=0;//initialize
  minnll=10e8;//initialize
  // while (stat>1){
	//   if (ntries>=3) break;
	  bern1_fit = bern1.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  // //   stat = bern1_fit->status();
	// //   minnll = bern1_fit->minNll();
	// //   if (stat>1) params_bern1->assignValueOnly(bern1_fit->randomizePars());
	// //   ntries++; 
	// // }
  // // stat=1;//initialize 
  // // ntries=0;//initialize
  // // minnll=10e8;//initialize
  // // while (stat>1){
	// //   if (ntries>=3) break;
	  bern2_fit = bern2.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  // //   stat = bern2_fit->status();
	// //   minnll = bern2_fit->minNll();
	// //   if (stat>1) params_bern2->assignValueOnly(bern2_fit->randomizePars());
	// //   ntries++; 
	// // }
  // // stat=1;//initialize 
  // // ntries=0;//initialize
  // // minnll=10e8;//initialize
  // // while (stat>1){
	// //   if (ntries>=3) break;
	  bern3_fit = bern3.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  // //   stat = bern3_fit->status();
	// //   minnll = bern3_fit->minNll();
	// //   if (stat>1) params_bern3->assignValueOnly(bern3_fit->randomizePars());
	// //   ntries++; 
	// // }
  // // stat=1;//initialize 
  // // ntries=0;//initialize
  // // minnll=10e8;//initialize
  // // while (stat>1){
	// //   if (ntries>=3) break;
	  bern4_fit = bern4.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  // //   stat = bern4_fit->status();
	// //   minnll = bern4_fit->minNll();
	// //   if (stat>1) params_bern4->assignValueOnly(bern4_fit->randomizePars());
	// //   ntries++; 
	// // }
  // // stat=1;//initialize 
  // // ntries=0;//initialize
  // // minnll=10e8;//initialize
  // // while (stat>1){
	// //   if (ntries>=3) break;
	  bern5_fit = bern5.fitTo(data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
  //   stat = bern5_fit->status();
	//   minnll = bern5_fit->minNll();
	//   if (stat>1) params_bern5->assignValueOnly(bern5_fit->randomizePars());
	//   ntries++; 
	// }
  

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
  
  RooWorkspace *ws =  new RooWorkspace();
  RooCategory catIndex(Form("pdfindex_ele_mu_cat%d_2020_13TeV",cat),Form("pdfindex_ele_mu_cat%d_2020_13TeV",cat));
  RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hzg_ele_mu_cat%d_2020_13TeV_bkgshape",cat),"all pdfs",catIndex,storedPdfs);
  RooRealVar nBackground(Form("CMS_hzg_ele_mu_cat%d_2020_13TeV_bkgshape_norm",cat),"nbkg",data.sumEntries(),0,10E8);
  ws->SetName("multipdf");
  ws->import(*pdf,RenameAllNodes(TString::Format("ele_mu_cat%d_2020", cat)), RenameAllVariablesExcept(TString::Format("ele_mu_cat%d_2020", cat), "mH"),RenameVariable("mH", "CMS_hzg_mass"));
  ws->import(nBackground,RenameVariable("BkgPdf_norm",Form("BkgPdf_ele_mu_cat%d_2020_norm",cat)));
  ws->import(catIndex);
  ws->import(data,RenameVariable("mH", "CMS_hzg_mass"));
  fout->cd();
  ws->Write();
  fout->Close();
  RooPlot* xframe1  = mH.frame() ;
  mH.setRange("blind1",100,120) ;
  mH.setRange("blind2",130,180);
  data.plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(80),CutRange("blind2")) ;
  gauxpow1.plotOn(xframe1,RooFit::Name("gauxpow1"),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashed));
  gauxpow3.plotOn(xframe1,RooFit::Name("gauxpow3"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashed));
  gauxpow5.plotOn(xframe1,RooFit::Name("gauxpow5"),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashed));
  gauxexp1.plotOn(xframe1,RooFit::Name("gauxexp1"),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDotted));
  gauxexp3.plotOn(xframe1,RooFit::Name("gauxexp3"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDotted));
  gauxexp5.plotOn(xframe1,RooFit::Name("gauxexp5"),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDotted));
  gauxlau1.plotOn(xframe1,RooFit::Name("gauxlau1"),LineColor(TColor::GetColor("#FF595E")),LineStyle(kDashDotted));
  gauxlau2.plotOn(xframe1,RooFit::Name("gauxlau2"),LineColor(TColor::GetColor("#1982C4")),LineStyle(kDashDotted));
  gauxlau3.plotOn(xframe1,RooFit::Name("gauxlau3"),LineColor(TColor::GetColor("#FFCA3A")),LineStyle(kDashDotted));
  bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(TColor::GetColor("#FF595E")));
  bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(TColor::GetColor("#FFCA3A")));
  bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(TColor::GetColor("#8AC926")));
  bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(TColor::GetColor("#1982C4")));
  bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(TColor::GetColor("#6A4C93")));

  
  // // int np_pow1 = gauxpow1.getParameters(data)->getSize();
  // // RooRealVar norm("norm","norm",data.sumEntries(),0,10E6);
  // // RooExtendPdf *pdf = new RooExtendPdf("ext","ext",gauxpow1,norm);
  // // double chi2_pow1 = xframe1->chiSquare("gauxpow1","data",np_pow1);
  // // cout<<"========> "<<endl;
  // // cout<<chi2_pow1<<endl;
  // // cout<<np_pow1<<endl;
  // // int np_pow3 = gauxpow1.getParameters(data)->getSize();
  // // double chi2_pow3 = xframe1->chiSquare("gauxpow3","data",np_pow3);
  // // cout<<chi2_pow3<<endl;
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
  leg4->AddEntry(xframe1->findObject("bern1"), "bern1", "l");
  leg4->AddEntry(xframe1->findObject("bern2"), "bern2", "l");
  leg4->AddEntry(xframe1->findObject("bern3"), "bern3", "l");
  leg4->AddEntry(xframe1->findObject("bern4"), "bern4", "l");
  leg4->AddEntry(xframe1->findObject("bern5"), "bern5", "l");
    leg4->AddEntry(xframe1->findObject("data"), "data", "lep");

  leg4->Draw("same");
  gPad->Print(Form("cat%d_turn.pdf",cat));
  cout<<"=======>Category"<<cat<<endl;
  if(lau1_fit->status()>1){cout<<"Lau1:"<<lau1_fit->status()<<endl;lau1_fit->Print();}
  if(lau2_fit->status()>1){cout<<"Lau2:"<<lau2_fit->status()<<endl;lau2_fit->Print();}
  if(lau3_fit->status()>1){cout<<"Lau3:"<<lau3_fit->status()<<endl;lau3_fit->Print();}
  if(pow1_fit->status()>1){cout<<"Pow1:"<<pow1_fit->status()<<endl;pow1_fit->Print();}
  if(pow3_fit->status()>1){cout<<"Pow3:"<<pow3_fit->status()<<endl;pow3_fit->Print();}
  if(pow5_fit->status()>1){cout<<"Pow5:"<<pow5_fit->status()<<endl;pow5_fit->Print();}
  if(exp1_fit->status()>1){cout<<"Exp1:"<<exp1_fit->status()<<endl;exp1_fit->Print();}
  if(exp3_fit->status()>1){cout<<"Exp3:"<<exp3_fit->status()<<endl;exp3_fit->Print();}
  if(exp5_fit->status()>1){cout<<"Exp5:"<<exp5_fit->status()<<endl;exp5_fit->Print();}
  if(bern1_fit->status()>1){cout<<"Bern1:"<<bern1_fit->status()<<endl;bern1_fit->Print();}
  if(bern2_fit->status()>1){cout<<"Bern2:"<<bern2_fit->status()<<endl;bern2_fit->Print();}
  if(bern3_fit->status()>1){cout<<"Bern3:"<<bern3_fit->status()<<endl;bern3_fit->Print();}
  if(bern4_fit->status()>1){cout<<"Bern4:"<<bern4_fit->status()<<endl;bern4_fit->Print();}
  if(bern5_fit->status()>1){cout<<"Bern5:"<<bern5_fit->status()<<endl;bern5_fit->Print();}
}
