#include "weightRemoveUnphysical2p2hExtendedEventsClass.h"

//using namespace PlotUtils;

PlotUtils::weightRemoveUnphysical2p2hExtendedEventsClass::weightRemoveUnphysical2p2hExtendedEventsClass(){
	// TODO constructor stub
}

PlotUtils::weightRemoveUnphysical2p2hExtendedEventsClass::~weightRemoveUnphysical2p2hExtendedEventsClass(){
	// TODO destructor stub
}

double PlotUtils::weightRemoveUnphysical2p2hExtendedEventsClass::get2p2h_mc_er_W(int mc_er_nPart,
					                              		 const std::vector<int>& mc_er_ID,
							                         const std::vector<int>& mc_er_status,
										 const std::vector<double>& mc_er_Px,
										 const std::vector<double>& mc_er_Py,
										 const std::vector<double>& mc_er_Pz,
										 const std::vector<double>& mc_er_E){

	
        double VmuVnu = 1.0;
	double Vmu = 1.0, Pmux = 1.0, Pmuy = 1.0, Pmuz = 1.0, Emu = 1.0;
	double Vnu = 1.0, Pnux = 1.0, Pnuy = 1.0, Pnuz = 1.0, Enu = 1.0;

        const double mev2gev = 0.001;
        const double MmuMeV = 105.6583715;
	const double MpMeV = 938.272;

	for(int i=0; i < mc_er_nPart; i++){
		if(mc_er_status[i]==0 && mc_er_ID[i]==14){
			Pnux = mc_er_Px[i]; Pnuy = mc_er_Py[i]; Pnuz = mc_er_Pz[i]; Enu = mc_er_E[i];
			Vnu=TMath::Sqrt(mc_er_Px[i]*mc_er_Px[i] + mc_er_Py[i]*mc_er_Py[i] + mc_er_Pz[i]*mc_er_Pz[i]);}
		if(mc_er_status[i]==1 && mc_er_ID[i]==13){
			Pmux = mc_er_Px[i]; Pmuy = mc_er_Py[i]; Pmuz = mc_er_Pz[i]; Emu = mc_er_E[i];
			Vmu=TMath::Sqrt(mc_er_Px[i]*mc_er_Px[i] + mc_er_Py[i]*mc_er_Py[i] + mc_er_Pz[i]*mc_er_Pz[i]);}

		VmuVnu = Pmux*Pnux + Pmuy*Pnuy + Pmuz*Pnuz;
	}
        
	double CosMu = VmuVnu/(Vmu*Vnu);
        
	double Q2 = 2.0*Enu*(Emu - Vmu*CosMu) - MmuMeV*MmuMeV;
	double W = TMath::Sqrt(MpMeV*MpMeV - Q2 + 2.0*MpMeV*(Enu-Emu));

	return W * mev2gev;

}

double PlotUtils::weightRemoveUnphysical2p2hExtendedEventsClass::getWeightInternal(double q0_t,
	                           						   int mc_intType,
			                              				   int mc_er_nPart,
					                              		   const std::vector<int>& mc_er_ID,
							                           const std::vector<int>& mc_er_status,
										   const std::vector<double>& mc_er_Px,
										   const std::vector<double>& mc_er_Py,
										   const std::vector<double>& mc_er_Pz,
										   const std::vector<double>& mc_er_E/*, int rw2p2h*/){


	//compute a rolloff of the unphysical portion of Nieves
	//do NOT apply this when weighting to SuSA, their prediction is ok.
	
	double W = get2p2h_mc_er_W(mc_er_nPart, mc_er_ID, mc_er_status, mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E);			
	double myweight = 1.0;

	if(/*rw2p2h == 5 &&*/ mc_intType == 8){
		if(W > 1.500){
			if(W > 1.700){
				myweight *= 0.0;
			}else {
				double dW1 = 1.700 - 1.500;
				double dW2 = 1.700 - W;
				myweight *= dW2 / dW1;
			}

			//double tq0 = q.Energy();
			//this also flattens the top a little bit.
			//if(q.Energy() > 1.800){
			if(q0_t > 1.800){
				if(q0_t < 2.000){
					double dq1 = 2.000 - 1.800;
					double dq2 = 2.000 - q0_t;
					myweight *= dq2 / dq1;
				}
			}
	      	}
	}
      	// optional, could transition more smoothly in q3
	//std::cout<<myweight<<std::endl;
      	return myweight;
 
  
}

double PlotUtils::weightRemoveUnphysical2p2hExtendedEventsClass::getWeight(double q0_t,
	                           					   int mc_intType,
			                              			   int mc_er_nPart,
					                              	   const std::vector<int>& mc_er_ID,
							                   const std::vector<int>& mc_er_status,
									   const std::vector<double>& mc_er_Px,
									   const std::vector<double>& mc_er_Py,
									   const std::vector<double>& mc_er_Pz,
									   const std::vector<double>& mc_er_E){
 
        return getWeightInternal(q0_t, mc_intType, mc_er_nPart, mc_er_ID, mc_er_status, mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E/*, int rw2p2h*/);			
}

