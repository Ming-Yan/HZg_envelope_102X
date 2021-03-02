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
void getmultipdf(int cat,vector<int> vec)
{
  TFile *f = TFile::Open(Form("turnon_cat%d.root",cat));
  RooWorkspace *w = (RooWorkspace*)f->Get("multipdf");
  RooRealVar *CMS_hzg_mass = (RooRealVar*) w->var("CMS_hzg_mass");
  RooDataSet *data = (RooDataSet*)w->data(Form("data_obs_ele_mu_cat%d_2020",cat));
  // RooCategory *cats = (RooCategory*)w->cat("catIndex");
  RooCategory catIndex("catIndex","c");
  RooRealVar nBackground("CMS_hzg_bkgshape_norm","nbkg",data->sumEntries(),0,10E8);
TString catname[13] = 
  {
    "gauxexp1", "gauxexp3", "gauxexp5",
    "gauxpow1", "gauxpow3", "gauxpow5",
    "gauxlau1", "gauxlau2", "gauxlau3",
    "bern2", "bern3", "bern4", "bern5",
  };
 RooArgList storedPdfs("store");
 TString color[5]={"#FF595E","#FFCA3A","#8AC926","#1982C4","#6A4C93"};
 TLegend* leg4 = new TLegend(0.6,0.6,0.9,0.9);
 leg4->SetTextSize(0.04);
leg4->SetNColumns(2);
RooPlot* xframe1  = CMS_hzg_mass->frame() ;
 for(int i = 0 ; i < vec.size(); i++)
{
  RooAbsPdf *plot = (RooAbsPdf*)w->pdf(catname[vec[i]].Data());
  storedPdfs.add(*plot);
  
  CMS_hzg_mass->setRange("blind1",100,120) ;
  CMS_hzg_mass->setRange("blind2",130,170);
  data->plotOn(xframe1,Binning(80),CutRange("blind1"),RooFit::Name("data"));
  data->plotOn(xframe1,Binning(80),CutRange("blind2"));
  int colorcode = 0;
  if(vec[i]==0||vec[i]==3||vec[i]==6) colorcode=0;
  else if(vec[i]==7||vec[i]==9) colorcode=1;
  else if(vec[i]==1||vec[i]==4||vec[i]==8||vec[i]==10)colorcode=2;
  else if(vec[i]==11) colorcode=3;
  else colorcode=4;
  if(catname[vec[i]].Contains("pow"))plot->plotOn(xframe1,RooFit::Name(catname[vec[i]].Data()),LineColor(TColor::GetColor(color[colorcode].Data())),LineStyle(kDashed));
  else if(catname[vec[i]].Contains("exp"))plot->plotOn(xframe1,RooFit::Name(catname[vec[i]].Data()),LineColor(TColor::GetColor(color[colorcode].Data())),LineStyle(kDotted));
  else if(catname[vec[i]].Contains("lau"))plot->plotOn(xframe1,RooFit::Name(catname[vec[i]].Data()),LineColor(TColor::GetColor(color[colorcode].Data())),LineStyle(kDashDotted));
  else if(catname[vec[i]].Contains("bern")) plot->plotOn(xframe1,RooFit::Name(catname[vec[i]].Data()),LineColor(TColor::GetColor(color[colorcode].Data())));
  leg4->AddEntry(xframe1->findObject(catname[vec[i]].Data()), catname[vec[i]].Data(), "l");
}
xframe1->SetMinimum(0.0001);
xframe1->Draw();
leg4->Draw("same");
gPad->Print(Form("turnon_ftest_cat%d.pdf",cat));
 RooMultiPdf *pdf = new RooMultiPdf("CMS_hzg_bkgshape","all pdfs",catIndex,storedPdfs);
 TFile *fout = new TFile(Form("new_cat%d.root",cat),"recreate");
 RooWorkspace *ws =  new RooWorkspace();
 ws->SetName("multipdf");
 ws->import(*pdf);
 ws->import(nBackground);
 ws->import(catIndex);
 ws->import(*data);
 fout->cd();
 ws->Write();
}
void multipdf()
{
  int cat[8]={1,2,3,4,501,502,503,6789};
  for(int i = 0 ; i < 8;i++){
    vector<int> vec;vec.clear();
    if(i==0){vec.push_back(0);vec.push_back(3);vec.push_back(6);vec.push_back(9);vec.push_back(10);}
    else if(i==1){vec.push_back(0);vec.push_back(3);vec.push_back(4);vec.push_back(6);vec.push_back(9);vec.push_back(10);}
    else if(i==2){vec.push_back(0);vec.push_back(3);vec.push_back(4);vec.push_back(6);vec.push_back(7);vec.push_back(8);vec.push_back(9);vec.push_back(10);}
    else if(i==3){vec.push_back(0);vec.push_back(1);vec.push_back(3);vec.push_back(6);vec.push_back(9);vec.push_back(10);}
    else if(i==4){vec.push_back(0);vec.push_back(3);vec.push_back(6);vec.push_back(9);}
    else if(i==5){vec.push_back(0);vec.push_back(3);vec.push_back(6);vec.push_back(7);vec.push_back(9);}
    else if(i==6){vec.push_back(0);vec.push_back(3);vec.push_back(6);vec.push_back(7);vec.push_back(9);vec.push_back(10);}
    else if(i==7){vec.push_back(0);vec.push_back(3);vec.push_back(6);vec.push_back(7);vec.push_back(9);}
    getmultipdf(cat[i],vec);
  }
  

}
