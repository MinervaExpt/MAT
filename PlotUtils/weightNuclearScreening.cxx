#include "weightNuclearScreening.h"
#include "TSystem.h"

using namespace PlotUtils;

void weightNuclearScreening::read(const std::string f)
{
  fNucScrFile = new TFile(f.c_str(),"READONLY");
  if(fNucScrFile) fNucScrRatio = (TH2D*)fNucScrFile->Get("nuclear");

}

double weightNuclearScreening::GetWeight(double nu,double Q2){//Inputs both in GeV!!
  
  if(nu>30) return 1; //Only have model up to 30 GeV
  if(Q2>1) return 1; //Only have model up to 1 GeV2

  int binnum = fNucScrRatio->FindBin(nu,Q2);
  return fNucScrRatio->GetBinContent(binnum);

}
