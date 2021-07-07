#include "weightMK.h"
#include "TreeWrapper.h"

using namespace PlotUtils;

//----------------------------------------------------------------------------------------------------
void weightMK::read(const TString  f)
//Read in the params doubles from a file
//The Weight is based in Single pion production in neutrino-nucleon interactions
//M. Kabirnezhad
//Phys. Rev. D 97, 013002 – Published 23 January 2018
//https://doi.org/10.1103/PhysRevD.97.013002

//argument: valid filename
{
  fMKratio = TFile::Open(f,"READONLY");
  if (fMKratio){

    hMKratio_CC1ppip = (TH2F*)fMKratio->Get("CC1ppip_rat");
    hMKratio_CC1ppip->Print();
    std::cout << "Have read CC1ppip ratios from file " << f <<std::endl;
    
    hMKratio_CC1pi0 = (TH2F*)fMKratio->Get("CC1pi0_rat");
    hMKratio_CC1pi0->Print();
    std::cout << "Have read CC1pi0 ratios from file " << f <<std::endl;
    
    hMKratio_CC1npip = (TH2F*)fMKratio->Get("CC1npip_rat");
    hMKratio_CC1npip->Print();
    std::cout << "Have read CC1npip ratios from file " << f <<std::endl;
    
  }
  else{
    std::cout << "File could not be read" << std::endl;
    
  }
}
//----------------------------------------------------------------------------------------------------
double weightMK::getWeightInternal(const double mc_W, const double mc_Q2, int variation){
  //assert(hMKratio_CC1ppip);
  Double_t thisrwtemp;

  if(mc_W < 1.7 && mc_W > 1.079 && mc_Q2 < 1.5 && mc_Q2 > 0.0){
	Int_t Wbin, Q2bin;

	switch (variation){
		case 0:	Wbin  = hMKratio_CC1ppip->GetXaxis()->FindBin(mc_W);
			Q2bin = hMKratio_CC1ppip->GetYaxis()->FindBin(mc_Q2);
			thisrwtemp = hMKratio_CC1ppip -> GetBinContent(Wbin,Q2bin);  
			break;
		case 1:	Wbin  = hMKratio_CC1pi0->GetXaxis()->FindBin(mc_W);
			Q2bin = hMKratio_CC1pi0->GetYaxis()->FindBin(mc_Q2);
			thisrwtemp = hMKratio_CC1pi0  -> GetBinContent(Wbin,Q2bin);  
			break;
		case 2:	Wbin  = hMKratio_CC1npip->GetXaxis()->FindBin(mc_W);
			Q2bin = hMKratio_CC1npip->GetYaxis()->FindBin(mc_Q2);
			thisrwtemp = hMKratio_CC1npip -> GetBinContent(Wbin,Q2bin);  
			break;
		default: thisrwtemp = 1.0;
			break;
		//default:  std::cout << "Unknown variation number "<< variation << ", returning 1.0. USE, 0: CC1ppip, 1: CC1pi0, 2: CC1npip."     << std::endl;
        }
  } 
  else thisrwtemp = 1.0; 
  
  return thisrwtemp;
}
//----------------------------------------------------------------------------------------------------

double weightMK::getWeight(const double mc_W, const double mc_Q2, int variation){
	return getWeightInternal(mc_W, mc_Q2, variation); 
}

//----------------------------------------------------------------------------------------------------
// Don't forget to use true GENIE values for mc_W anf mc_Q2
// and they also should be in GeV

double weightMK::getWeight(const double mc_W,  // should true GENIE mc_w and in GeV
	                   const double mc_Q2, // should true GENIE mc_Q2 and in GeV 
			   const int mc_current, 
			   const int mc_intType, 
                           const int mc_er_nPart, 
			   const int mc_targetA, 
			   const int *mc_er_ID, 
			   const int *mc_er_status)
{
	
	int channel=-1, pdg, status, nn=0, np=0, npi0=0, npip=0, npISt=0, nnISt=0;
        bool isCCRES = mc_current == 1 && (mc_intType == 2 || (mc_intType==3 && mc_W < 1.7));

	if(isCCRES){
		if(mc_targetA == 1){
			for(int i=0; i<mc_er_nPart; i++){
				pdg = mc_er_ID[i];
				status = mc_er_status[i];
				
				if(status==0 && pdg==2212) npISt++;

				if(status==1 && pdg==2212) np++;     
				if(status==1 && pdg==2112) nn++;     
				if(status==1 && pdg==211 ) npip++;     
				if(status==1 && pdg==111 ) npi0++;	
			}
		}
		else{
			for(int i=0; i<mc_er_nPart; i++){
				pdg = mc_er_ID[i];
				status = mc_er_status[i];
				
				if(status==11 && pdg==2212) npISt++;
				if(status==11 && pdg==2112) nnISt++;

				if(status==14 && pdg==2212) np++;     
				if(status==14 && pdg==2112) nn++;     
				if(status==14 && pdg==211 ) npip++;     
				if(status==14 && pdg==111 ) npi0++;	
			}
		}
		
		if(npISt==1 && nnISt==0 && np==1 && nn==0 && npip==1 && npi0==0) channel=0;
		if(npISt==0 && nnISt==1 && np==1 && nn==0 && npip==0 && npi0==1) channel=1;
		if(npISt==0 && nnISt==1 && np==0 && nn==1 && npip==1 && npi0==0) channel=2;

		return getWeight(mc_W, mc_Q2, channel);
	}
	else return 1.0;
}

//----------------------------------------------------------------------------------------------------
// Don't forget to use true GENIE values for mc_W anf mc_Q2
// and they also should be in GeV

double weightMK::getWeight(const double mc_W,  // should true GENIE mc_w and in GeV
	                   const double mc_Q2, // should true GENIE mc_Q2 and in GeV 
			   const int mc_current, 
			   const int mc_intType, 
                           const int mc_er_nPart, 
			   const int mc_targetA, 
			   const std::vector<int>& mc_er_ID, 
			   const std::vector<int>& mc_er_status)
{
	
	int channel=-1, pdg, status, nn=0, np=0, npi0=0, npip=0, npISt=0, nnISt=0;
        bool isCCRES = mc_current == 1 && (mc_intType == 2 || (mc_intType==3 && mc_W < 1.7));

	if(isCCRES){
		if(mc_targetA == 1){
			for(int i=0; i<mc_er_nPart; i++){
				pdg = mc_er_ID[i];
				status = mc_er_status[i];
				
				if(status==0 && pdg==2212) npISt++;

				if(status==1 && pdg==2212) np++;     
				if(status==1 && pdg==2112) nn++;     
				if(status==1 && pdg==211 ) npip++;     
				if(status==1 && pdg==111 ) npi0++;	
			}
		}
		else{
			for(int i=0; i<mc_er_nPart; i++){
				pdg = mc_er_ID[i];
				status = mc_er_status[i];
				
				if(status==11 && pdg==2212) npISt++;
				if(status==11 && pdg==2112) nnISt++;

				if(status==14 && pdg==2212) np++;     
				if(status==14 && pdg==2112) nn++;     
				if(status==14 && pdg==211 ) npip++;     
				if(status==14 && pdg==111 ) npi0++;	
			}
		}
		
		if(npISt==1 && nnISt==0 && np==1 && nn==0 && npip==1 && npi0==0) channel=0;
		if(npISt==0 && nnISt==1 && np==1 && nn==0 && npip==0 && npi0==1) channel=1;
		if(npISt==0 && nnISt==1 && np==0 && nn==1 && npip==1 && npi0==0) channel=2;

		return getWeight(mc_W, mc_Q2, channel);
	}
	else return 1.0;
}
//----------------------------------------------------------------------------------------------------
double weightMK::getWeight(PlotUtils::TreeWrapper* chw,
			   Long64_t entry)
{

	const double mev2gev  = 0.001;
	const double mc_Q2    = mev2gev*mev2gev*chw->GetValue("mc_Q2", entry);
	const double mc_w     = mev2gev*chw->GetValue("mc_w", entry);
	const int mc_current  = chw->GetValue("mc_current", entry); 
	const int mc_intType  = chw->GetValue("mc_intType", entry); 
	const int mc_targetA  = chw->GetValue("mc_targetA", entry); 
	const int mc_er_nPart = chw->GetValue("mc_er_nPart",entry);

	std::vector<int> mc_er_ID     = chw->GetValueVector<int>("mc_er_ID", entry); 
	std::vector<int> mc_er_status = chw->GetValueVector<int>("mc_er_status", entry); 	

	double ret = 1.0;

        ret = getWeight(mc_w,
		        mc_Q2,
			mc_current,
			mc_intType,
			mc_er_nPart,
			mc_targetA,
			mc_er_ID,
			mc_er_status);

	return ret;
}


PlotUtils::weightMK& PlotUtils::weight_mk() {
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string dir_data = std::string(mparalocation)+"/data/Reweight/";
  static PlotUtils::weightMK* _weight_MK = 
      new PlotUtils::weightMK(dir_data+"output_ratio_genie_neut_for_MKmodel.root");
  return *_weight_MK;
}    
