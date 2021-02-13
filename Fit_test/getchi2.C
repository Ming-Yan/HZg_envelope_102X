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
void getchi2(int cat)
{
  TFile *f = TFile::Open(Form("turnon_cat%d.root",cat));
  RooWorkspace *w = (RooWorkspace*)f->Get("multipdf");
    string catname[14] = 
      {
        "bern1", "bern2", "bern3", "bern4", "bern5",
        "gauxexp1", "gauxexp3", "gauxexp5",
        "gauxlau1", "gauxlau2", "gauxlau3",
        "gauxpow1", "gauxpow3", "gauxpow5"
      };
    
    RooRealVar *CMS_hzg_mass = (RooRealVar*) w->var("CMS_hzg_mass");
    // CMS_hzg_mass.setBins(320);
    RooPlot* xframe1  = CMS_hzg_mass->frame() ;
    RooDataSet *data = (RooDataSet*)w->data(Form("data_obs_ele_mu_cat%d_2020",cat));
    vector <RooAbsPdf*> pdf; pdf.clear();
    ofstream outtxt(Form("chi2_cat%d.txt",cat));
    cout<<"====> "<<cat<<" <===="<<endl;
    outtxt<<"Functions\tndf\tchi2/ndf"<<endl;
    for(int i = 0; i < 14; i++)
      {
        RooAbsPdf *pdf = (RooAbsPdf*)w->pdf(catname[i].c_str());
        // pdf->Print();
        data->plotOn(xframe1,Binning(320),RooFit::Name("data")) ;
        int np = pdf->getParameters(data)->getSize();
        RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
        pdf->plotOn(xframe1,RooFit::Name(catname[i].c_str()));
        RooExtendPdf *pdf_ext = new RooExtendPdf("ext","ext",*pdf,norm);
        double chi2 = xframe1->chiSquare(catname[i].c_str(),"data",np);
        cout<<catname[i].c_str()<<" "<<chi2<<endl;
        outtxt<<catname[i].c_str()<<"\t"<<np<<"\t"<<chi2<<endl;
      }
}
