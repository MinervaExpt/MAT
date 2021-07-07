#ifndef weightRemoveUnphysical2p2hExtendedEventsClass_h
#define weightRemoveUnphysical2p2hExtendedEventsClass_h

#include <iostream> //cout
#include <TMath.h>
#include <vector>
#include <Rtypes.h>

namespace PlotUtils{

	class weightRemoveUnphysical2p2hExtendedEventsClass{
	public:
	  weightRemoveUnphysical2p2hExtendedEventsClass();

	  virtual ~weightRemoveUnphysical2p2hExtendedEventsClass(); 

	  double getWeight(double q0_t,
			   int mc_intType,
			   int mc_er_nPart,
			   const std::vector<int>& mc_er_ID,
			   const std::vector<int>& mc_er_status,
			   const std::vector<double>& mc_er_Px,
			   const std::vector<double>& mc_er_Py,
			   const std::vector<double>& mc_er_Pz,
			   const std::vector<double>& mc_er_E);

	  private:
	  double get2p2h_mc_er_W(int mc_er_nPart,
				  const std::vector<int>& mc_er_ID,
				  const std::vector<int>& mc_er_status,
				  const std::vector<double>& mc_er_Px,
				  const std::vector<double>& mc_er_Py,
				  const std::vector<double>& mc_er_Pz,
				  const std::vector<double>& mc_er_E);
	  
	  double getWeightInternal(double q0_t,
			  	   int mc_intType,
			  	   int mc_er_nPart,
				   const std::vector<int>& mc_er_ID,
				   const std::vector<int>& mc_er_status,
				   const std::vector<double>& mc_er_Px,
				   const std::vector<double>& mc_er_Py,
				   const std::vector<double>& mc_er_Pz,
				   const std::vector<double>& mc_er_E);
	  
	};
}
#endif
