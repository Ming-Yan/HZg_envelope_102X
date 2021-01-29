#ifndef H_SETCONSTPARAM
#define H_SETCONSTPARAM
using namespace RooFit;
void SetConstantParams(const RooArgSet* params)
{
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next())
    {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      if (rrv)
        {
	  rrv->setConstant(true);
	  // std::cout << " " << rrv->GetName();
        }
    }
} // close set const parameters
#endif
