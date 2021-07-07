#include "weight_minervaq2qe.h"

using namespace PlotUtils;

void weight_minervaq2qe::read(TString filename)
//Read in the params doubles from a file
//argument: valid filename
{

  fQ2QERatio = TFile::Open(filename,"READONLY");
  if(fQ2QERatio){
    hQ2QERatio = (TH1D*)fQ2QERatio->Get("mc_q2qe_CV_WithStatErr");
  }
  else{
    std::cout << "Bad file input for weight_minervaq2qe. Try again" << std::endl;
    exit(1);
  }
}




double weight_minervaq2qe::getWeightInternal(const double q2qe) {

  double retval = 1;
  double checkval = q2qe;
  if(q2qe>10) checkval=9;
  if(q2qe<0){
    std::cout << "You have a q2qe passed to weight_minervaq2qe less than 0. Non-physical. Returning 1.0"<< std::endl;
    return 1.0;
  }
  retval = hQ2QERatio->Interpolate(checkval);
  return retval;
}


double weight_minervaq2qe::getWeight(const double q2qe) {
  return getWeightInternal(q2qe);
}






