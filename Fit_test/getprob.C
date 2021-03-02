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

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  //params_test->Print("v");
  int stat=1;
  double minnll=10e8;
  while (stat!=0){
    if (ntries>=MaxTries) break;
    RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1)
				       ,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;
}
int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){


  double global_minNll = 1E10;
  int best_index = 0;
  int number_of_indeces = cat->numTypes();
		
  RooArgSet snap,clean;
  RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
  params->remove(*cat);
  params->snapshot(snap);
  params->snapshot(clean);
  if (!silent) {
    params->Print("V");
  }
	
  //bkg->setDirtyInhibit(1);
  //RooAbsReal *nllm = bkg->createNLL(*data);
  //RooMinimizer minim(*nllm);
  //minim.setStrategy(1);
	
  for (int id=0;id<number_of_indeces;id++){		
    params->assignValueOnly(clean);
    cat->setIndex(id);

    //RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

    if (!silent) {
      
	std::cout << "BEFORE  MAKING FIT" << std::endl;
	params->Print("V");
	std::cout << "-----------------------" << std::endl;		
	
    }
		
    //minim.minimize("Minuit2","minimize");
    double minNll=0; //(nllm->getVal())+bkg->getCorrection();
    int fitStatus=1;		
    runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/5);
  
    minNll=minNll+bkg->getCorrection();

    if (!silent) {
      
	std::cout << "After Minimization ------------------  " <<std::endl;
	std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
	bkg->Print("v");
	bkg->getCurrentPdf()->getParameters(*data)->Print("V");
	std::cout << " ------------------------------------  " << std::endl;
	
	params->Print("V");
      
      std::cout << "[INFO] AFTER FITTING" << std::endl;
      std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
      std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
      std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
      std::cout << "-----------------------" << std::endl;
    }
			
    if (minNll < global_minNll){
      global_minNll = minNll;
      snap.assignValueOnly(*params);
      best_index=id;
    }
  }
  cat->setIndex(best_index);
  params->assignValueOnly(snap);
	
  if (!silent) {
  std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
  //bkg->getCurrentPdf()->getParameters(*data)->Print("v");
  }
  return best_index;
}

void getprob(int cat)
{
  TFile *f = TFile::Open(Form("turnon_cat%d.root",cat));
  RooWorkspace *w = (RooWorkspace*)f->Get("multipdf");
    string catname[13] = 
      {
        "gauxexp1", "gauxexp3", "gauxexp5",
        "gauxpow1", "gauxpow3", "gauxpow5",
        "gauxlau1", "gauxlau2", "gauxlau3",
        "bern2", "bern3", "bern4", "bern5",
      };
    RooRealVar *CMS_hzg_mass = (RooRealVar*) w->var("CMS_hzg_mass");
    RooCategory *cats = (RooCategory*)w->cat("catIndex");
    // CMS_hzg_mass.setBins(320);
    RooPlot* xframe1  = CMS_hzg_mass->frame() ;
    RooDataSet *data = (RooDataSet*)w->data(Form("data_obs_ele_mu_cat%d_2020",cat));
    RooMultiPdf *pdf = (RooMultiPdf*)w->pdf("CMS_hzg_bkgshape");
     cout<<"Get best fit functions"<<endl;
      int bestfit = 0;
      bestfit=getBestFitFunction(pdf, data, cats);
      cout<<catname[bestfit]<<endl;
    for(int i = 0 ; i < 13; i++)
      {
        RooFitResult *result = (RooFitResult*)f->Get(Form("fitresult_%s_data_obs_ele_mu_cat%d_2020",catname[i].c_str(),cat));
        double prevNLL,thisNLL,chi2,prob;
        int order,prevorder;
        if(i<3){
          order = i*2+1;
          thisNLL = result->minNll();
          chi2 = 2.*(prevNLL-thisNLL);
          if (chi2<0. && order>1) chi2=0.;
          if(i>0)
          {
            prob = TMath::Prob(chi2,order-prevorder);
            cout<<catname[i]<<" "<<prob<<endl;
          }
          if(i<2){prevNLL = thisNLL; prevorder=order;}
        }//Exponential
        else if(i<6){
          order = (i-3)*2+1;
          thisNLL = result->minNll();
          chi2 = 2.*(prevNLL-thisNLL);
          if (chi2<0. && order>1) chi2=0.;
          if(i>3)
          {
            prob = TMath::Prob(chi2,order-prevorder);
            cout<<catname[i]<<" "<<prob<<endl;
          }
          if(i<5){prevNLL = thisNLL; prevorder=order;}
        }//PowerLaw
        else if(i<9){
          order=i-5;
          thisNLL = result->minNll();
          chi2 = 2.*(prevNLL-thisNLL);
          if (chi2<0. && order>1) chi2=0.;
          if(i>6)
          {
            prob = TMath::Prob(chi2,order-prevorder);
            cout<<catname[i]<<" "<<prob<<endl;
          }
          if(i<8){prevNLL = thisNLL; prevorder=order;}        
          }//Laurent
        else {
          order=i-7;
          thisNLL = result->minNll();
          chi2 = 2.*(prevNLL-thisNLL);
          if (chi2<0. && order>1) chi2=0.;
          if(i>9)
          {
            prob = TMath::Prob(chi2,order-prevorder);
            cout<<catname[i]<<" "<<prob<<endl;
          }
          if(i<12){prevNLL = thisNLL; prevorder=order;}
        }//Bernstein
        
      }
     
}
