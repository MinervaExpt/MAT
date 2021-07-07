#ifndef weightGenieBodekRitchieClass_h
#define weightGenieBodekRitchieClass_h

#include <iostream> //cout
#include <TMath.h>
#include <vector>
#include <Rtypes.h>

namespace PlotUtils{

	class weightGenieBodekRitchieClass{
	public:
	  weightGenieBodekRitchieClass();
	  virtual ~weightGenieBodekRitchieClass(); 
	  double getWeight(int rwBRtail,
			   int mc_er_nPart,
			   int mc_intType,
			   int mc_targetA,
			   const std::vector<int>& mc_er_status,
			   const std::vector<int>& mc_er_ID,
			   const std::vector<double>& mc_er_Px,
			   const std::vector<double>& mc_er_Py,
			   const std::vector<double>& mc_er_Pz,
			   bool verbose = false); //in GeV

	  private:
	  double getWeightInternal(int rwBRtail,
				   int mc_er_nPart,
				   int mc_intType,
				   int mc_targetA,
				   const std::vector<int>& mc_er_status,
				   const std::vector<int>& mc_er_ID,
				   const std::vector<double>& mc_er_Px,
				   const std::vector<double>& mc_er_Py,
				   const std::vector<double>& mc_er_Pz,
				   bool verbose);
	  
	  
	};
}
#endif
