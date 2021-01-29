#include "RooStepBernstein.h"
#include "RooGaussStepBernstein.h"
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
  RooDecay exp("exp","", mH,tau,turnon,RooDecay::SingleSided);
  RooRealVar meanp("meanp","meanp",0);
  RooRealVar sigmap("sigmap","sigmap",2,0.01,10);
  RooRealVar alphap("alphap","alphap",105,50,200);
  RooRealVar betap("betap","betap",6,0,20);
  RooGaussModel turnonp("turnonp","",mH,meanp,sigmap);
  RooGenericPdf tail ("tail","tail","1e-20 + (@0 > @1)*((@0)^(-@2))",RooArgList(mH,alphap,betap));
  RooFFTConvPdf pow("pow","pow", mH, tail, turnonp);


  //RooGaussStepBernstein bern6("bern6","bern6",mH,mean,sigma,step, RooArgList(p0,p1,p2,p3,p4,p5,p6));
  
  /*RooRealVar edge("edge"," ",115, 100,120);
  RooGenericPdf stepFuncPDF("stepFuncPDF", "stepFuncPDF", "(@0=> @1)", RooArgList(mH,edge));
  
  RooRealVar cexp("cexp","cexp",-0.01) ;
  RooExponential exp("exp","exp",mH,cexp);
  RooProdPdf prod2("prod2", "prod2",RooArgList(exp,stepFuncPDF));
  RooRealVar sigma("sigma", "Resolution Model Sigma",8.0);
  RooRealVar mean("mean", "mean",100);
  RooGaussian resMod("resMod", "Resolution Model", mH, mean, sigma);
  mH.setBins(40000, "cache");

  RooFFTConvPdf sxg("sxg", "step (X) gauss", mH,prod2,resMod);
  */  
    RooFitResult* bern2_fit= bern2.fitTo(data, Range("window"), Save(kTRUE)) ;
    RooFitResult* bern1_fit= bern1.fitTo(data, Range("window"), Save(kTRUE)) ;

  RooFitResult* bern3_fit= bern3.fitTo(data, Range("window"), Save(kTRUE)) ;
  RooFitResult* bern4_fit= bern4.fitTo(data, Range("window"), Save(kTRUE)) ;
  RooFitResult* bern5_fit= bern5.fitTo(data, Range("window"), Save(kTRUE)) ;
  RooFitResult* exp_fit = exp.fitTo(data, Range("window"), Save(kTRUE));

  RooFitResult* pow_fit =pow.fitTo(data, Range("window"), Save(kTRUE));
  //RooFitResult* bern6_fit= bern6.fitTo(data, Range("window"), Save(kTRUE)) ;
  RooPlot* xframe1  = mH.frame() ;
  mH.setRange("blind1",100,120) ;
  mH.setRange("blind2",130,180);
  data.plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(80),CutRange("blind2")) ;
  //prod.plotOn(xframe1);
    bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(kGreen));

    bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(kMagenta));

  bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(kBlue));
  bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(kRed));
  bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(kOrange-3));
  exp.plotOn(xframe1,RooFit::Name("exp"),LineColor(kAzure+1),LineStyle(kDashed));
  pow.plotOn(xframe1,RooFit::Name("pow"),LineColor(kGray+2),LineStyle(kDashed));
  //bern6.plotOn(xframe1, RooFit::Name("bern6"),LineColor(kGreen-2));
  //stepFuncPDF.plotOn(xframe1);

  //sxg.plotOn(xframe1,RooFit::Name("sxg"),LineColor(kRed));        
  xframe1->SetMinimum(0.0001);
  xframe1->Draw();
  TLegend* leg4 = new TLegend(0.7,0.5,0.9,0.9);
  leg4->AddEntry(xframe1->findObject("data"), "data", "lep");
  leg4->AddEntry(xframe1->findObject("bern1"), "bern1", "l");
  leg4->AddEntry(xframe1->findObject("bern2"), "bern2", "l");
  leg4->AddEntry(xframe1->findObject("bern3"), "bern3", "l");
  leg4->AddEntry(xframe1->findObject("bern4"), "bern4", "l");
  leg4->AddEntry(xframe1->findObject("bern5"), "bern5", "l");
  leg4->AddEntry(xframe1->findObject("exp"), "exp", "l");
  leg4->AddEntry(xframe1->findObject("pow"), "pow", "l");
  //leg4->AddEntry(xframe1->findObject("bern6"), "bern6", "l");
  leg4->Draw("same");

}
