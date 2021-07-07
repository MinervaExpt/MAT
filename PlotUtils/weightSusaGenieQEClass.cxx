#include "weightSusaGenieQEClass.h"

using namespace PlotUtils;

void weightSusaGenieQEClass::read(const TString  f)
//Read in the params doubles from a file
//argument: valid filename
{
  fSusaGenieQE = TFile::Open(f,"READONLY");
  if (fSusaGenieQE){

    hRatio = (TH2D*)fSusaGenieQE->Get("hRatio");
    
    if(NULL == hRatio)std::cout << "hRatioPN is null " << std::endl;
    hRatio->Print();
    std::cout << "have read in ratios from file " << f <<std::endl;
    
  } else{
    //File could not be read
    std::cout << "File could not be read" << std::endl;
    
  }
}

 
double weightSusaGenieQEClass::getWeightInternal(const double mc_q0,
	       				         const double mc_q3,
						 bool verbose){

  //std::cout << "assert hRatio" << std::endl;
  if(NULL == hRatio)std::cout << "hRatio is NULL" << std::endl;
  assert(hRatio);
  
  // the construction here is that the histogram bins
  // line up in mev-sized steps e.g. from 0.018 to 0.019
  // and the value stored is the value at bin-center.
  
  int q3bin = hRatio->GetXaxis()->FindBin(mc_q3);
  int q0bin = hRatio->GetYaxis()->FindBin(mc_q0);

  // There is a 27 MeV offset baked into the ratio.
  // There are no events below 0.017, but just in case.
  double thisrwtemp = 1.0, weight=1.0;
  if(mc_q3 < 2.0){
	  thisrwtemp = hRatio->GetBinContent(q3bin,q0bin);

	  if(verbose)std::cout << "SuSA QE weight " << thisrwtemp << " q0q3 " << mc_q0 << " " << mc_q3 << " " << " q0q3 " << q0bin << " " << q3bin << std::endl;
	  if(thisrwtemp >= 0.0 && thisrwtemp < 100.0) weight *= thisrwtemp;
	  else {
		std::cout << " Unphysical Weight? " << thisrwtemp << " q0q3 " << mc_q0 << " " << mc_q3 << " q0q3 " << q0bin << " " << q3bin  << std::endl;
		weight = 1.0;}

	  return weight;
  } else return 1.0;
}

double weightSusaGenieQEClass::getWeightMeV(const double mc_q0, const double mc_q3, bool verbose){

  return getWeightInternal(mc_q0*0.001, mc_q3*0.001, verbose);

}

double weightSusaGenieQEClass::getWeight(const double mc_q0,
	       			         const double mc_q3,
					 bool verbose){

  return getWeightInternal(mc_q0, mc_q3, verbose);

}



