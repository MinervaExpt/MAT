#include "weightDIS.h"

using namespace PlotUtils;

void weightDIS::read(enum modelSet m)
//Read in the params doubles from a file
//argument: valid filename
{
  char *mparalocation = std::getenv("MPARAMFILESROOT");

  //assumes modelSet==knCTEQ15
  std::string f_C = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQ15_C_total_sigma_GENIE_C.root";
  std::string f_Fe = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQ15_Fe_total_sigma_GENIE_Fe.root";
  std::string f_Pb = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQ15_Pb_total_sigma_GENIE_Pb.root";
  if(m==weightDIS::knCTEQ15){
    std::cout << "I will be loading the nCTQE15 reweighting values" << std::endl;
  }
  else if(m==weightDIS::kAMU){
    std::cout << "I will be loading the AMU reweighting values" << std::endl;
    f_C = std::string(mparalocation)+"/data/Reweight/AMU_C_GENIE.root";;
    f_Fe = std::string(mparalocation)+"/data/Reweight/AMU_Fe_GENIE.root";;
    f_Pb = std::string(mparalocation)+"/data/Reweight/AMU_Pb_GENIE.root";;
  }
  else if(m==weightDIS::knCTEQ15nu){
    f_C = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQnu_C_total_sigma_GENIE_C.root";
    f_Fe = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQnu_Fe_total_sigma_GENIE_Fe.root";
    f_Pb = std::string(mparalocation)+"/data/Reweight/total_ratios_total_sigma_nCTEQnu_Pb_total_sigma_GENIE_Pb.root";
  }
  else{
    std::cout << "You didn't specify a known modelSet. "
                 "Please consult Ana/MCReweight/MCReweight/weightDIS.h for valid modelSets. "
                 "If you modelset isn't list, but you think it should contact software experts." 
              << std::endl;
    exit(1);
  }
  


  fDISRatio_C = TFile::Open(f_C.c_str(),"READONLY");
  fDISRatio_Fe = TFile::Open(f_Fe.c_str(),"READONLY");
  fDISRatio_Pb = TFile::Open(f_Pb.c_str(),"READONLY");
  if (fDISRatio_C){
    hDISRatio_C = (TH3D*)fDISRatio_C->Get("h_weights");
    std::cout << "have read in ratios from file " << f_C <<std::endl;
  }
  else{
    //File could not be read
    std::cout << "C file could not be read" << std::endl;
  }

  if (fDISRatio_Fe){
    hDISRatio_Fe = (TH3D*)fDISRatio_Fe->Get("h_weights");
    std::cout << "have read in ratios from file " << f_Fe <<std::endl;
  }
  else{
    //File could not be read
    std::cout << "Fe file could not be read" << std::endl;
  }

  if (fDISRatio_Pb){
    hDISRatio_Pb = (TH3D*)fDISRatio_Pb->Get("h_weights");
    std::cout << "have read in ratios from file " << f_Pb <<std::endl;
  }
  else{
    //File could not be read
    std::cout << "Pb file could not be read" << std::endl;
  }
}


double weightDIS::getWeightInternal(const double x, const double y, 
                                    const double enu, const int targetZ) {
  assert(hDISRatio_C);
  assert(hDISRatio_Fe);
  assert(hDISRatio_Pb);
  double retval = 1;

  //Figure out what bin in xyz we are getting
  int binx = hDISRatio_C->GetXaxis()->FindBin(x); //xbin is Bjorken x
  int biny = hDISRatio_C->GetYaxis()->FindBin(enu); //ybin is neutrino energy
  int binz = hDISRatio_C->GetZaxis()->FindBin(y); //zbin is inelasticity

  if(targetZ==82)      retval = hDISRatio_Pb->GetBinContent(binx,biny,binz);//Lead
  else if(targetZ==26) retval = hDISRatio_Fe->GetBinContent(binx,biny,binz);//Iron
  else if(targetZ==1)  retval = 1.0;//hydrogen has no weighting
  else                 retval = hDISRatio_C->GetBinContent(binx,biny,binz);//Carbon plus everyting else remaining (Oxygen mostly)

  if(retval < 1e-3 || retval > 100){
    if (!m_quiet){
      std::cout << "weightDIS has an extreme event weight with the following values (x,y,Enu,targetZ)" 
		<< x << "\t" << y << "\t" << enu << "\t" << targetZ
		<< ". The weight value is " << retval << std::endl;
      std::cout << "Setting weight to 1." << std::endl;
    }
    retval=1.0;
  }
  return retval;
}


double weightDIS::getWeight(const double x, const double y, const double enu,
                            const int targetZ) {
  return getWeightInternal(x,y,enu,targetZ);
}






