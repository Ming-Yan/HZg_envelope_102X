using namespace RooFit;
void getchi2()
{
  TFile *f = TFile::Open("higgsCombineacceff_1_bf.MultiDimFit.mH125.38.root");
  RooWorkspace *w  = (RooWorkspace *)f->Get("w");
  RooDataSet *data = (RooDataSet*)w->data("data_obs");
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  RooRealVar *mass = w->var("CMS_hzg_mass");
  RooPlot *plot_chi2 = mass->frame();

  //RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);
  //w.pdf("model_s"), w.pdf("model_b")
  RooAbsPdf *pdf = w->pdf("model_b");
  //RooAbsPdf *pdf = pdfa->getPdf("cat1");
  //pdf->Print();
  data->plotOn(plot_chi2,Binning(220),Name("data"));
  pdf->plotOn(plot_chi2,Name("cat1"));
  int np = pdf->getParameters(*data)->getSize();
  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  double prob = TMath::Prob(chi2*(220-np),220-np);
  cout<<prob<<endl;
}
