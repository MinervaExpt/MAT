#include "weightGenieBodekRitchieClass.h"

//using namespace PlotUtils;

PlotUtils::weightGenieBodekRitchieClass::weightGenieBodekRitchieClass(){
	// TODO constructor stub
}

PlotUtils::weightGenieBodekRitchieClass::~weightGenieBodekRitchieClass(){
	// TODO destructor stub
}
double PlotUtils::weightGenieBodekRitchieClass::getWeightInternal(	int rwBRtail,
								int mc_er_nPart,
								int mc_intType,
								int mc_targetA,
								const std::vector<int>& mc_er_status,
				                		const std::vector<int>& mc_er_ID,
								const std::vector<double>& mc_er_Px,
								const std::vector<double>& mc_er_Py,
								const std::vector<double>& mc_er_Pz,
								bool verbose){
 
  double gevmev = 0.001 ;
  double weight=1.0; 
  double struckNucleonPdg, struckNucleonP; 

  for(int i=0; i < mc_er_nPart; ++i){
  	if(mc_er_status[i] == 11){
		struckNucleonPdg = mc_er_ID[i];
		struckNucleonP   = TMath::Sqrt(mc_er_Px[i]*mc_er_Px[i] + mc_er_Py[i]*mc_er_Py[i] + mc_er_Pz[i]*mc_er_Pz[i]);
	}
  }

  double weightBRtail = 1.0;

  if(mc_intType==1 && mc_targetA==12 && (struckNucleonPdg==2212||struckNucleonPdg==2112)){
  if(verbose)std::cout <<"struckNucleon"<< struckNucleonPdg << " " << struckNucleonP << " " << weightBRtail << " " << mc_targetA << " " << mc_intType << std::endl;
  	// Pb and Fe and Ti sometimes have their mc_er arrays corrupted
  	// Can't do those weights anyway.
  	// But the low threshold is nucleus dependent, so come back to this.

  	double lowthreshold  = 0.221/gevmev;
  	double highthreshold = 0.500/gevmev;
  
	if(struckNucleonP > lowthreshold){
  		// at threshold weight up by 6.   at 500 no additional weight.
  		//double factor = 5.0;
  		double factor = 5.0;
  		double diff =  1.0 - ( struckNucleonP - lowthreshold)/(highthreshold-lowthreshold);
  		//weightBRtail = 1.0 + factor*diff;
  		weightBRtail = 1.0 + factor*diff;
	        if(verbose)std::cout <<"lowthreshold "<< lowthreshold <<" highthreshold "<< highthreshold <<" weightBRtail "<< weightBRtail << std::endl;
  	} else {
  	// mode 2, also weight down the main part.
  		if(rwBRtail==2)weightBRtail = 0.74;
  		else weightBRtail = 1.0;
  	}
	
  	if(verbose)std::cout << "struckNucleon " << struckNucleonPdg << " " << struckNucleonP << " " << weightBRtail << std::endl;
  }
  weight *= weightBRtail;
  
  return weight;
  
}

double PlotUtils::weightGenieBodekRitchieClass::getWeight( int rwBRtail,
		     			int mc_er_nPart,
					int mc_intType,
					int mc_targetA,
					const std::vector<int>& mc_er_status,
			                const std::vector<int>& mc_er_ID,
					const std::vector<double>& mc_er_Px,
					const std::vector<double>& mc_er_Py,
					const std::vector<double>& mc_er_Pz,
					bool verbose){

  return getWeightInternal(rwBRtail, mc_er_nPart, mc_intType, mc_targetA, mc_er_status, mc_er_ID, mc_er_Px, mc_er_Py, mc_er_Pz, verbose);

}



