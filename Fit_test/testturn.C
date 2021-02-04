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
#include "/data3/mingyan/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
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
  RooArgList storedPdfs("store");
  TFile *fout  = new TFile("turnon_leptag.root","recreate");
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
  RooRealVar mean_pow3("mean_pow3","mean_pow3",0);
  RooRealVar sigma_pow3("sigma_pow3","sigma_pow3",2,0.0001,20);
  RooRealVar turnon_pow3("turnon_pow3","turnon_pow3",105,80,120);
  RooRealVar mean_pow5("mean_pow5","mean_pow5",0);
  RooRealVar sigma_pow5("sigma_pow5","sigma_pow5",2,0.0001,20);
  RooRealVar turnon_pow5("turnon_pow5","turnon_pow5",105,80,120);
  RooRealVar p1_pow1("p1_pow1","p1_pow1",0.5,-10,10.);
  RooRealVar cp1_pow1("cp1_pow1","cp1_pow1",0.5,0,1.);
  RooRealVar p1_pow3("p1_pow3","p1_pow3",0.5,-10,10.);
  RooRealVar cp1_pow3("cp1_pow3","cp1_pow3",1.,0.,1.);
  RooRealVar p3_pow3("p3_pow3","p3_pow3",0.5,-10.,10.);
  RooRealVar cp3_pow3("cp3_pow3","cp3_pow3",0.5,0,1.);
  RooRealVar p1_pow5("p1_pow5","p1_pow5",0.5,-10,10.);
  RooRealVar cp1_pow5("cp1_pow5","cp1_pow5",1.,0.,1.);
  RooRealVar p3_pow5("p3_pow5","p3_pow5",0.5,-10.,10.);
  RooRealVar cp3_pow5("cp3_pow5","cp3_pow5",0.5,0,1.);
  RooRealVar p5_pow5("p5_pow5","p5_pow5",0.5,-10.,10.);
  RooRealVar cp5_pow5("cp5_pow5","cp5_pow5",0.5,0,1.);
  RooGenericPdf step_pow1("step_pow1", "step_pow1", "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(mH,turnon_pow1,p1_pow1,cp1_pow1));//step*(ax^b)
  RooGenericPdf step_pow3("step_pow3", "step_pow3", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(mH,turnon_pow3,p1_pow3,cp1_pow3,p3_pow3,cp3_pow3));//step*(ax^b+cx^d)
  RooGenericPdf step_pow5("step_pow5", "step_pow5", "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(mH,turnon_pow5,p1_pow5,cp1_pow5,p3_pow5,cp3_pow5,p5_pow5,cp5_pow5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_pow1("gau_pow1","gau_pow1",mH,mean_pow1,sigma_pow1);
  RooGaussModel gau_pow3("gau_pow3","gau_pow3",mH,mean_pow3,sigma_pow3);
  RooGaussModel gau_pow5("gau_pow5","gau_pow5",mH,mean_pow5,sigma_pow5);
  RooFFTConvPdf gauxpow1("gauxpow1","gauxpow1",mH,step_pow1,gau_pow1);
  RooFFTConvPdf gauxpow3("gauxpow3","gauxpow3",mH,step_pow3,gau_pow3);
  RooFFTConvPdf gauxpow5("gauxpow5","gauxpow5",mH,step_pow5,gau_pow5);
  
  //testing with generic Laurent
  RooRealVar mean_lau1("mean_lau1","mean_lau1",0);
  RooRealVar sigma_lau1("sigma_lau1","sigma_lau1",2,0.0001,20);
  RooRealVar turnon_lau1("turnon_lau1","turnon_lau1",105,80,120);
  RooRealVar mean_lau2("mean_lau2","mean_lau2",0);
  RooRealVar sigma_lau2("sigma_lau2","sigma_lau2",2,0.0001,20);
  RooRealVar turnon_lau2("turnon_lau2","turnon_lau2",105,80,120);
  RooRealVar mean_lau3("mean_lau3","mean_lau3",0);
  RooRealVar sigma_lau3("sigma_lau3","sigma_lau3",2,0.0001,20);
  RooRealVar turnon_lau3("turnon_lau3","turnon_lau3",105,80,120);
  RooRealVar cl1_lau1("cl1_lau1","cl1_lau1",0.25,0,1.);
  RooRealVar cl2_lau1("cl2_lau1","cl2_lau1",0.25,0,1.);
  RooRealVar cl1_lau2("cl1_lau2","cl1_lau2",0.25,0.,1.);
  RooRealVar cl2_lau2("cl2_lau2","cl2_lau2",0.25,0.,1.);
  RooRealVar cl3_lau2("cl3_lau2","cl3_lau2",0.25/2.,0,1.);
  RooRealVar cl1_lau3("cl1_lau3","cl1_lau3",0.25,0.,1.);
  RooRealVar cl2_lau3("cl2_lau3","cl2_lau3",0.25,0.,1.);
  RooRealVar cl3_lau3("cl3_lau3","cl3_lau3",0.25/2.,0,1.);
  RooRealVar cl4_lau3("cl4_lau3","cl4_lau3",0.25/3.,0,1.);
  RooGenericPdf step_lau1("step_lau1", "step_lau1", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(mH,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
  RooGenericPdf step_lau2("step_lau2", "step_lau2", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(mH,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
  RooGenericPdf step_lau3("step_lau3", "step_lau3", "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(mH,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
  RooGaussModel gau_lau1("gau_lau1","gau_lau1",mH,mean_lau1,sigma_lau1);
  RooGaussModel gau_lau2("gau_lau2","gau_lau2",mH,mean_lau2,sigma_lau2);
  RooGaussModel gau_lau3("gau_lau3","gau_lau3",mH,mean_lau3,sigma_lau3);
  RooFFTConvPdf gauxlau1("gauxlau1","gauxlau1",mH,step_lau1,gau_lau1);
  RooFFTConvPdf gauxlau2("gauxlau2","gauxlau2",mH,step_lau2,gau_lau2);
  RooFFTConvPdf gauxlau3("gauxlau3","gauxlau3",mH,step_lau3,gau_lau3);


  // testing with generic exponential
  RooRealVar mean_exp1("mean_exp1","mean_exp1",0);
  RooRealVar sigma_exp1("sigma_exp1","sigma_exp1",2,0.0001,20);
  RooRealVar turnon_exp1("turnon_exp1","turnon_exp1",105,80,120);
  RooRealVar mean_exp3("mean_exp3","mean_exp3",0);
  RooRealVar sigma_exp3("sigma_exp3","sigma_exp3",2,0.0001,20);
  RooRealVar turnon_exp3("turnon_exp3","turnon_exp3",105,80,120);
  RooRealVar mean_exp5("mean_exp5","mean_exp5",0);
  RooRealVar sigma_exp5("sigma_exp5","sigma_exp5",2,0.0001,20);
  RooRealVar turnon_exp5("turnon_exp5","turnon_exp5",105,80,120);
  RooRealVar p1_exp1("p1_exp1","p1_exp1",-0.01,-0.5,0.);
  RooRealVar cp1_exp1("cp1_exp1","cp1_exp1",0.5,0,1.);
  RooRealVar p1_exp3("p1_exp3","p1_exp3",-0.01,-0.5,0.);
  RooRealVar cp1_exp3("cp1_exp3","cp1_exp3",1.,0.,1.);
  RooRealVar p3_exp3("p3_exp3","p3_exp3",-0.01,-0.5,0.);
  RooRealVar cp3_exp3("cp3_exp3","cp3_exp3",0.5,0,1.);
  RooRealVar p1_exp5("p1_exp5","p1_exp5",-0.01,-0.5,0.);
  RooRealVar cp1_exp5("cp1_exp5","cp1_exp5",1.,0.,1.);
  RooRealVar p3_exp5("p3_exp5","p3_exp5",-0.01,-0.5,0.);
  RooRealVar cp3_exp5("cp3_exp5","cp3_exp5",0.5,0,1.);
  RooRealVar p5_exp5("p5_exp5","p5_exp5",-0.01,-0.5,0.);
  RooRealVar cp5_exp5("cp5_exp5","cp5_exp5",0.5,0,1.);
  RooGenericPdf step_exp1("step_exp1", "step_exp1", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(mH,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
  RooGenericPdf step_exp3("step_exp3", "step_exp3", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(mH,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
  RooGenericPdf step_exp5("step_exp5", "step_exp5", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(mH,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
  RooGaussModel gau_exp1("gau_exp1","gau_exp1",mH,mean_exp1,sigma_exp1);
  RooGaussModel gau_exp3("gau_exp3","gau_exp3",mH,mean_exp3,sigma_exp3);
  RooGaussModel gau_exp5("gau_exp5","gau_exp5",mH,mean_exp5,sigma_exp5);
  RooFFTConvPdf gauxexp1("gauxexp1","gauxexp1",mH,step_exp1,gau_exp1);
  RooFFTConvPdf gauxexp3("gauxexp3","gauxexp3",mH,step_exp3,gau_exp3);
  RooFFTConvPdf gauxexp5("gauxexp5","gauxexp5",mH,step_exp5,gau_exp5);

  RooFitResult *pow1_fit = gauxpow1.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *pow3_fit = gauxpow3.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *pow5_fit = gauxpow5.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *exp1_fit = gauxexp1.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *exp3_fit = gauxexp3.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *exp5_fit = gauxexp5.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *lau1_fit = gauxlau1.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *lau2_fit = gauxlau2.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *lau3_fit = gauxlau3.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *bern1_fit = bern1.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *bern2_fit = bern2.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *bern3_fit = bern3.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *bern4_fit = bern4.fitTo(data,Range("window"),Save(kTRUE));
  RooFitResult *bern5_fit = bern5.fitTo(data,Range("window"),Save(kTRUE));

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
  RooCategory catIndex("catIndex","c");
  RooMultiPdf *pdf = new RooMultiPdf("CMS_hzg_bkgshape","all pdfs",catIndex,storedPdfs);
  RooRealVar nBackground("CMS_hzg_bkgshape_norm","nbkg",data.sumEntries(),0,10E8);
  ws->SetName("multipdf");
  ws->import(*pdf);
  ws->import(nBackground);
  ws->import(catIndex);
  ws->import(data);
  fout->cd();
  ws->Write();
  fout->Close();
  RooPlot* xframe1  = mH.frame() ;
  mH.setRange("blind1",100,120) ;
  mH.setRange("blind2",130,180);
  data.plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data")) ;
  data.plotOn(xframe1,Binning(80),CutRange("blind2")) ;
  // gauxpow1.plotOn(xframe1,RooFit::Name("gauxpow1"));
  // gauxpow3.plotOn(xframe1,RooFit::Name("gauxpow3"),LineColor(kOrange-3));
  // gauxpow5.plotOn(xframe1,RooFit::Name("gauxpow5"),LineColor(kRed));
  gauxexp1.plotOn(xframe1,RooFit::Name("gauxexp1"));
  gauxexp3.plotOn(xframe1,RooFit::Name("gauxexp3"),LineColor(kOrange-3));
  gauxexp5.plotOn(xframe1,RooFit::Name("gauxexp5"),LineColor(kRed));
  // gauxlau1.plotOn(xframe1,RooFit::Name("gauxlau1"));
  // gauxlau2.plotOn(xframe1,RooFit::Name("gauxlau2"),LineColor(kRed));
  // gauxlau3.plotOn(xframe1,RooFit::Name("gauxlau3"),LineColor(kOrange-3));
  // bern1.plotOn(xframe1, RooFit::Name("bern1"),LineColor(kGreen));
  // bern2.plotOn(xframe1, RooFit::Name("bern2"),LineColor(kMagenta));
  // bern3.plotOn(xframe1, RooFit::Name("bern3"),LineColor(kBlue));
  // bern4.plotOn(xframe1, RooFit::Name("bern4"),LineColor(kRed));
  // bern5.plotOn(xframe1, RooFit::Name("bern5"),LineColor(kOrange-3));
  // exp.plotOn(xframe1,RooFit::Name("exp"),LineColor(kAzure+1),LineStyle(kDashed));
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
  // leg4->AddEntry(xframe1->findObject("gauxpow1"), "pow1", "l");
  // leg4->AddEntry(xframe1->findObject("gauxpow3"), "pow3", "l");
  // leg4->AddEntry(xframe1->findObject("gauxpow5"), "pow5", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp1"), "exp1", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp3"), "exp3", "l");
  leg4->AddEntry(xframe1->findObject("gauxexp5"), "exp5", "l");
  // leg4->AddEntry(xframe1->findObject("gauxlau1"), "lau1", "l");
  // leg4->AddEntry(xframe1->findObject("gauxlau2"), "lau2", "l");
  // leg4->AddEntry(xframe1->findObject("gauxlau3"), "lau3", "l");
  // //leg4->AddEntry(xframe1->findObject("bern6"), "bern6", "l");
  leg4->Draw("same");

}
