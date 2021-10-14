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
void getchi2(int cat,int md, int mu)
{
  //int mu=170;int md=105;
  TFile *f = TFile::Open(Form("m%d_%d_turnon_cat%d_new.root",md,mu,cat));
  //TFile *f = TFile::Open(Form("rename_turnon_cat%d.root",cat));
  RooWorkspace *w = (RooWorkspace*)f->Get("multipdf");
  //w->Print();
    string catname[16] = 
      {
	//"gauxexp1","gauxpow1","gauxlau3","bern2","bern3"
	"bern1", "bern2", "bern3", "bern4", "bern5",
        "gauxexp1", "gauxexp3", "gauxexp5",
        "gauxlau1", "gauxlau2", "gauxlau3","gauxlau4","gauxlau5",
        "gauxpow1", "gauxpow3", "gauxpow5"
      };
    //w->Print(); 
    RooRealVar *CMS_hzg_mass = (RooRealVar*) w->var("CMS_hzg_mass");
    // CMS_hzg_mass.setBins(320);


    vector <RooAbsPdf*> pdf; pdf.clear();
    ofstream outtxt(Form("chi2_cat%d.txt",cat));
    cout<<"====> "<<cat<<" <===="<<endl;
    outtxt<<"Functions\tndf\tchi2/ndf"<<endl;
    RooPlot* xframe1  = CMS_hzg_mass->frame() ;
    //RooDataHist *data = (RooDataHist*)w->data(Form("CMS_hzg_datahist_ele_mu_cat%d_2020_13TeV",cat));
    //RooDataSet *data = (RooDataSet*)w->data("data_obs_ele_mu_cat4_2020");
    RooDataSet *data = (RooDataSet*)w->data(Form("data_obs_ele_mu_cat%d_2020_13TeV",cat));
    for(int i =  1; i <16; i++)
      {
	RooAbsPdf *pdf = (RooAbsPdf*)w->pdf(Form("%s_env_pdf_ele_mu_cat%d_2020_13TeV",catname[i].c_str(),cat));
	//RooAbsPdf *pdf = (RooAbsPdf*)w->pdf(catname[i].c_str());
	///pdf->Print();
        data->plotOn(xframe1,Binning(4*(mu-md)),RooFit::Name("data")) ;
        int np = pdf->getParameters(data)->getSize();
        RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
        pdf->plotOn(xframe1,RooFit::Name(catname[i].c_str()));
        RooExtendPdf *pdf_ext = new RooExtendPdf("ext","ext",*pdf,norm);
	//RooFitResult *result = (RooFitResult*)f->Get(Form("fitresult_%s_datahist_ele_mu_cat%d_2020",catname[i].c_str(),cat));
        double chi2 = xframe1->chiSquare(catname[i].c_str(),"data",np);
	//pdf->plotOn(xframe1,VisualizeError(*result,1),FillColor(kOrange)) ;
	//xframe1->Draw();
	//gPad->Print(Form("cat%d_%s.C",cat,catname[i].c_str()));
        cout<<catname[i].c_str()<<" "<<np<<" "<<chi2<<" "<<TMath::Prob(chi2*(4*(mu-md)-np),4*(mu-md)-np)<<endl;
	//outtxt<<catname[i].c_str()<<"\t"<<np<<"\t"<<chi2<<" "<<TMath::Prob(chi2,np)<<endl;
	}
}
