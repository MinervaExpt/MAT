#include "weightSusaValenciaClass.h"

using namespace PlotUtils;

void weightSusaValenciaClass::read(const TString  f)
//Read in the params doubles from a file
//argument: valid filename
{
  fSusaValencia = TFile::Open(f,"READONLY");
  if (fSusaValencia){

    //hRatio = (TH2D*)fSusaValencia->Get("hRatioSusaValencia");
    hRatio = (TH2D*)fSusaValencia->Get("hRatio");
    //hRatioPN = (TH2D*)fSusaValencia->Get("hRatioSusaValenicaPN");
    hRatioPN = (TH2D*)fSusaValencia->Get("hRatioPN");
    //hRatioNotPN = (TH2D*)fSusaValencia->Get("hRatioSusaValenciaNotPN");
    hRatioNotPN = (TH2D*)fSusaValencia->Get("hRatioNotPN");


    if(NULL == hRatioPN)std::cout << "hRatioPN is null " << std::endl;
    hRatioPN->Print();
    std::cout << "have read in ratios from file " << f <<std::endl;
    
  }
  else{
    //File could not be read
    std::cout << "File could not be read" << std::endl;
    
  }
}

 
double weightSusaValenciaClass::getWeightInternal(const double mc_q0, const double mc_q3, const int mode) const{
  //assert(hSusaValenciaPN);
  
  // the construction here is that the histogram bins
  // line up in mev-sized steps e.g. from 0.018 to 0.019
  // and the value stored is the value at bin-center.
  // Mode 
  // 0 : Ratio
  // 1 : RatioPN
  // 2 : RatioNotPN

  int q3bin, q0bin;
  double thisrwtemp = 1.0;
  if(mc_q3 < 2.0){
	  switch(mode){
		  case 0 : q3bin = hRatio->GetXaxis()->FindBin(mc_q3);
			   q0bin = hRatio->GetYaxis()->FindBin(mc_q0);
			   thisrwtemp = hRatio->GetBinContent(q3bin,q0bin);
			   break;
		  case 1 : q3bin = hRatioPN->GetXaxis()->FindBin(mc_q3);
			   q0bin = hRatioPN->GetYaxis()->FindBin(mc_q0);
			   thisrwtemp = hRatioPN->GetBinContent(q3bin,q0bin);
			   break;
		  case 2 : q3bin = hRatioNotPN->GetXaxis()->FindBin(mc_q3);
			   q0bin = hRatioNotPN->GetYaxis()->FindBin(mc_q0);
			   thisrwtemp = hRatioNotPN->GetBinContent(q3bin,q0bin);
			   break;
		  default : thisrwtemp = 1.0;
			   break;
	  }

	  if(thisrwtemp >= 0.0 && thisrwtemp < 100.0) return thisrwtemp;
	  std::cout << " Unphysical Weight? " << thisrwtemp << std::endl;
  	  return 0.0;
  }
  else { return thisrwtemp;}
  //return 0.0; 
}

double weightSusaValenciaClass::getWeightMeV(const double mc_q0, const double mc_q3, const int mode) const{

  return getWeightInternal(mc_q0*0.001, mc_q3*0.001, mode);

}

double weightSusaValenciaClass::getWeight(const double mc_q0, const double mc_q3, const int mode) const{

  return getWeightInternal(mc_q0, mc_q3, mode);

}



