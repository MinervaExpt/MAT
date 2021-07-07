#ifndef weightSusaGenieQEClass_h
#define weightSusaGenieQEClass_h
//#endif

/*
memos to myself
sanitize the input histogram, in addition to the test here.
why did long titles in the histogram names not work?
make naming conventions exact between file (put SusaValenciaPN and NotPN in  the names)
change the capitalization?
add the high W rolloff here too, as an extra function and option.
add both GeV and MeV options
*/

#include <fstream>  //ifstream
#include <iostream> //cout

#include <TString.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "math.h"
#include "assert.h"
//For Compatibility with ROOT compiler
//uncomment the following:
//#define ROOT

/*!
 *  Code example to extract the RPA effect central weight
 *  and its uncertainties from the prepared files.
 *  Heidi Schellman (Oregon State) and Rik Gran (Minnesota Duluth)
 *  for use in MINERvA experiment analysis
 *  must compile with the ROOT libraries
 *  g++ `root-config --glibs --cflags` -O3 weightRPAtest.cxx -o weightRPAtest
 

 *    The underlying model is from the IFIC Valencia group
 *    see (public) minerva docdb:12688 for the full physics discussion

 */


// NOTE UNITS ARE GEV in the calculation
// make sure you convert MeV to GeV before calling these functions


// Class for getting  RPA paramaters inputs from a given file
// Includes methods that return all five RPA weights at once
//   (is the most cpu efficient way to get them)
// Or return just the central value
//   (skipping the uncertainties code completely if only cv wanted)
// Or return each CV and uncertainty one at a time
//   (is the least cpu efficient, repeats some calculations 5 times)
namespace PlotUtils{
	class weightSusaGenieQEClass {
	public:
	  //Constructor: Read in params from a filename
	  weightSusaGenieQEClass(const TString  f) { read(f); }    //Read in params from file
	  
	  TString  filename;
	  TFile* fSusaGenieQE;
	  TH2D *hRatio;
	  TH2D *hRatioPN;
	  TH2D *hRatioNotPN;
	  
	  //TH2D *hRPArelratio;
	  //TH2D *hRPAnonrelratio;
	  //TH1D *hQ2relratio;
	  //TH1D *hQ2nonrelratio;
	  //TArrayD *TADrpapolyrel;
	  //Double_t *rpapolyrel ;
	  //TArrayD *TADrpapolynonrel;
	  //Double_t *rpapolynonrel;
	  //static const int CENTRAL=0;
	  //static const int LOWQ2 = 1;
	  //static const int HIGHQ2 = 2;

	  // true to take from histogram, false use parameterization
	  //static const bool Q2histORparam = true;  

	  // MINERvA holds kinematics in MeV, but all these functions require GeV
	  // So make sure you pass them in GeV.
	  double getWeightMeV(const double mc_q0,
			      const double mc_q3,
			      bool verbose); //in GeV

	  double getWeight(const double mc_q0,
			   const double mc_q3,
			   bool verbose = false); //in GeV
	  //Initializer
	  void read(const TString  f);
	  
	  // q0 and q3 in GeV, type = 1 for low Q2, 2 for high Q2, 0 for central
	  //double getWeightInternal(const double mc_q0, const double mc_q3,int type, int sign);

	  private:
	  //double getWeightInternal(const double mc_Q2);
	  double getWeightInternal(const double mc_q0,
				   const double mc_q3,
				   bool verbose = false);
	  
	  
	};
}
#endif

