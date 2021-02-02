#include "RooStepBernstein.h"
#include "RooGaussStepBernstein.h"
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
void testturn()
{
  TTree* DataTree1 = makeTTree("HiggsMass_ptwei_bkg_ele_mu_cat6789_2020_full.txt");
  int totalev = DataTree1->GetEntriesFast();
  RooRealVar mH("mH", "mH", 100, 180, "GeV") ;
  mH.setRange("window",100, 180);
  RooDataSet data("data", " ", RooArgSet(mH), Import(*DataTree1));
  //generetic Bernstein polynomials
  RooRealVar mean("mean","mean",0);
  RooRealVar sigma("sigma","sigma",3,0.01,20);
  RooRealVar step("step","step",110,100,120);
  RooRealVar p0("p0","p0",15);
  RooRealVar b1p1("b1p1","b1p1",0.3,-1e-6,900);
  RooRealVar b2p1("b2p1","b2p1",0.3,-1e-6,900);
  RooRealVar b3p1("b3p1","b3p1",0.3,-1e-6,900);
  RooRealVar b4p1("b4p1","b4p1",0.3,-1e-6,900);
  RooRealVar b5p1("b5p1","b5p1",0.3,-1e-6,900);
  RooRealVar b2p2("b2p2","b2p2",0.3,-1e-6,900);
  RooRealVar b3p2("b3p2","b3p2",0.3,-1e-6,900);
  RooRealVar b4p2("b4p2","b4p2",0.3,-1e-6,900);
  RooRealVar b5p2("b5p2","b5p2",0.3,-1e-6,900);
  RooRealVar b3p3("b3p3","b3p3",0.3,-1e-6,900);
  RooRealVar b4p3("b4p3","b4p3",0.3,-1e-6,900);
  RooRealVar b5p3("b5p3","b5p3",0.3,-1e-6,900);
  RooRealVar b3p4("b3p4","b3p4",0.3,-1e-6,900);
  RooRealVar b4p4("b4p4","b4p4",0.3,-1e-6,900);
  RooRealVar b5p4("b5p4","b5p4",0.3,-1e-6,900);
  RooRealVar b3p5("b3p5","b3p5",0.3,-1e-6,900);
  RooRealVar b4p5("b4p5","b4p5",0.3,-1e-6,900);
  RooRealVar b5p5("b5p5","b5p5",0.3,-1e-6,900);
  RooGaussStepBernstein bern1("bern1","bern1",mH,mean,sigma,step, RooArgList(p0,b1p1));
  RooGaussStepBernstein bern2("bern2","bern2",mH,mean,sigma,step, RooArgList(p0,b2p1,b2p2));
  RooGaussStepBernstein bern3("bern3","bern3",mH,mean,sigma,step, RooArgList(p0,b3p1,b3p2,b3p3));
  RooGaussStepBernstein bern4("bern4","bern4",mH,mean,sigma,step, RooArgList(p0,b4p1,b4p2,b4p3,b4p4));
  RooGaussStepBernstein bern5("bern5","bern5",mH,mean,sigma,step, RooArgList(p0,b5p1,b5p2,b5p3,b5p4,b5p5));
  RooRealVar meang("meang","meang",120,90,150);
  RooRealVar sigmag("sigmag","sigmag",1,0.01,10);
  RooRealVar tau("tau","tau",5,0,50);
  RooGaussModel turnon("turnon","",mH,meang,sigmag);
    //testing with generic power law 
    RooRealVar mean_pow1("mean_pow1","mean_pow1",0);
    RooRealVar sigma_pow1("sigma_pow1","sigma_pow1",2,0.0001,20);
    RooRealVar turnon_pow1("turnon_pow1","turnon_pow1",105,80,120);
    RooRealVar mean_pow2("mean_pow2","mean_pow2",0);
    RooRealVar sigma_pow2("sigma_pow2","sigma_pow2",2,0.0001,20);
    RooRealVar turnon_pow2("turnon_pow2","turnon_pow2",105,80,120);
    RooRealVar mean_pow3("mean_pow3","mean_pow3",0);
    RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",2,0.0001,20);
    RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",105,80,120);
    RooRealVar p1_pow1("p1_pow1","p1_pow1",0.5,-10,10.);
    RooRealVar p1_pow2("p1_pow2","p1_pow2",0.5,-10,10.);
    RooRealVar p1_pow3("p1_pow3","p1_pow3",0.5,-10,10.);
    RooRealVar p2_pow2("p2_pow2","p2_pow2",1.,-10,10.);
    RooRealVar p2_pow3("p2_pow3","p2_pow3",1.,-10,10.);
    RooRealVar p3_pow3("p3_pow3","p3_pow3",0.2,-10,10.);
    RooGenericPdf step_pow1("step_pow1", "step_pow1", "1e-20+(@0> @1)*((@2/(@0*@0)))", RooArgList(mH,turnon_pow1,p1_pow1));//step*(ax^-2)
    RooGenericPdf step_pow2("step_pow2", "step_pow2", "1e-20+(@0> @1)*((@2/(@0*@0))+(@3/(@0*@0*@0)))", RooArgList(mH,turnon_pow2,p1_pow2,p2_pow2));//step*(ax^-2+bx^-3)
    RooGenericPdf step_pow3("step_pow3", "step_pow3", "1e-20+(@0> @1)*((@2/(@0*@0))+(@3/(@0*@0*@0))+(@4/(@0*@0*@0*@0)))", RooArgList(mH,turnon_pow3,p1_pow3,p2_pow3,p3_pow3));//step*(ax^-2+bx^-3+cx^-4)
    RooGaussModel gau_pow1("gau_pow1","gau_pow1",mH,mean_pow1,sigma_pow1);
    RooGaussModel gau_pow2("gau_pow2","gau_pow2",mH,mean_pow2,sigma_pow2);
    RooGaussModel gau_pow3("gau_pow3","gau_pow3",mH,mean_pow3,sigma_pow3);
    RooFFTConvPdf gauxpow1("gauxpow1","gauxpow1",mH,step_pow1,gau_pow1);
    RooFFTConvPdf gauxpow2("gauxpow2","gauxpow2",mH,step_pow2,gau_pow2);
    RooFFTConvPdf gauxpow3("gauxpow3","gauxpow3",mH,step_pow3,gau_pow3);

  // RooExponential expo("expo", "",mH)
  // RooDecay exp("exp","", mH,tau,turnon,RooDecay::SingleSided);
    //   RooFitResult* bern2_fit= bern2.fitTo(data, Range("window"), Save(kTRUE)) ;
  //   RooFitResult* bern1_fit= bern1.fitTo(data, Range("window"), Save(kTRUE)) ;

  // RooFitResult* bern3_fit= bern3.fitTo(data, Range("window"), Save(kTRUE)) ;
  // RooFitResult* bern4_fit= bern4.fitTo(data, Range("window"), Save(kTRUE)) ;
  // RooFitResult* bern5_fit= bern5.fitTo(data, Range("window"), Save(kTRUE)) ;
  // RooFitResult* exp_fit = exp.fitTo(data, Range("window"), Save(kTRUE));

  // RooFitResult* pow_fit =pow->fitTo(data, Range("window"), Save(kTRUE));
  //RooFitResult* bern6_fit= bern6.fitTo(data, Range("window"), Save(kTRUE)) ;
  RooFitResult *pow1_fit = gauxpow1.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *pow2_fit = gauxpow2.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *pow3_fit = gauxpow3.fitTo(data,Range("window"),Save(kTRUE));
  RooPlot* xframe1  = mH.frame() ;
  mH.setRange("blind1",100,120) ;
  mH.setRange("blind2",130,180);
  // gauxpow.plotOn(xframe1);
  // stepxpow.plotOn(xframe1);
  // pow->plotOn(xframe1);
 
  // step_pow.plotOn(xframe1);
  // tail.plotOn(xframe1,LineColor(kRed));
  data.plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(80),CutRange("blind2")) ;
  // //prod.plotOn(xframe1);
  gauxpow1.plotOn(xframe1,RooFit::Name("gauxpow1"));
  gauxpow2.plotOn(xframe1,RooFit::Name("gauxpow2"),LineColor(kRed));
  gauxpow3.plotOn(xframe1,RooFit::Name("gauxpow3"),LineColor(kOrange-3));
  //   bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(kGreen));
  //   bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(kMagenta));
  // bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(kBlue));
  // bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(kRed));
  // bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(kOrange-3));
  // exp.plotOn(xframe1,RooFit::Name("exp"),LineColor(kAzure+1),LineStyle(kDashed));
  // pow->plotOn(xframe1,RooFit::Name("pow"),LineColor(kGray+2),LineStyle(kDashed));
  // powgau.plotOn(xframe1,LineColor(kBlue));    
  xframe1->SetMinimum(0.0001);
  xframe1->Draw();
  TLegend* leg4 = new TLegend(0.7,0.5,0.9,0.9);
  leg4->AddEntry(xframe1->findObject("data"), "data", "lep");
  // leg4->AddEntry(xframe1->findObject("bern1"), "bern1", "l");
  // leg4->AddEntry(xframe1->findObject("bern2"), "bern2", "l");
  // leg4->AddEntry(xframe1->findObject("bern3"), "bern3", "l");
  // leg4->AddEntry(xframe1->findObject("bern4"), "bern4", "l");
  // leg4->AddEntry(xframe1->findObject("bern5"), "bern5", "l");
  // leg4->AddEntry(xframe1->findObject("exp"), "exp", "l");
  leg4->AddEntry(xframe1->findObject("gauxpow1"), "pow1", "l");
  leg4->AddEntry(xframe1->findObject("gauxpow2"), "pow2", "l");
  leg4->AddEntry(xframe1->findObject("gauxpow3"), "pow3", "l");
  // //leg4->AddEntry(xframe1->findObject("bern6"), "bern6", "l");
  leg4->Draw("same");

}