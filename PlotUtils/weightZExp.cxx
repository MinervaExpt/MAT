#include "weightZExp.h"

using namespace PlotUtils;

void weightZExp::read(const TString  f)
//Read in the params doubles from a file
//argument: valid filename
{
  fZExpRatioCov = TFile::Open(f,"READONLY");
  if (fZExpRatioCov){
    hZExpRatioCov = (PlotUtils::MnvH1D*)fZExpRatioCov->Get("zexpvar");
    std::cout << "have read in ratios from file " << f <<std::endl;
    hZExpRatioBand = hZExpRatioCov->GetVertErrorBand("zexpCov");
  }
  else{
    //File could not be read
    std::cout << "File could not be read" << std::endl;
  }
}


double weightZExp::getWeightInternal(const double mc_Q2) {
  assert(hZExpRatioCov);
  
  if(mc_Q2>8) return hZExpRatioCov->Interpolate(8);
  else return hZExpRatioCov->Interpolate(mc_Q2);
}


double weightZExp::getWeight(const double mc_Q2) {
  return getWeightInternal(mc_Q2);
}


double weightZExp::getWeightInternal(const double mc_Q2, int uni_n) {
  assert(hZExpRatioBand);
  TH1D *universe = (TH1D*)hZExpRatioBand->GetHist(uni_n);

  if(mc_Q2>8) return universe->Interpolate(8);
  else return universe->Interpolate(mc_Q2);

}

double weightZExp::getWeight(const double mc_Q2, int uni_n) {
  return getWeightInternal(mc_Q2, uni_n);
}




