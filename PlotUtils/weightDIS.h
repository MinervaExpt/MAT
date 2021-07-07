#ifndef weightDIS_h
#define weightDIS_h

#include <fstream>  //ifstream
#include <iostream> //cout

#include <TString.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "math.h"
#include "assert.h"
#include "MnvH1D.h"
#include "MnvVertErrorBand.h"
//For Compatibility with ROOT compiler
//uncomment the following:
//#define ROOT

/*!
  This is a class for taking a x,y,Enu value, plus a target number to reweight the GENIE DIS model to either nCTEQ15 or the AMU model.
TODO: Add reference to nCTEQ15 and AMU model here
 */

namespace PlotUtils{

  class weightDIS
  {
    public:
      enum modelSet{
        knCTEQ15,//Add ref
        kAMU,//Add ref
	knCTEQ15nu//Add ref
      };
      
      bool m_quiet;
      
      //Constructor: Read in params from a filename
      weightDIS(enum modelSet m, bool quiet=false) { read(m); m_quiet=quiet;}    //Read in params from files
      
      TString  filename;
      TFile* fDISRatio_C;
      TFile* fDISRatio_Fe;
      TFile* fDISRatio_Pb;
      TH3D *hDISRatio_C;
      TH3D *hDISRatio_Fe;
      TH3D *hDISRatio_Pb;

      // MINERvA holds kinematics in MeV, but all these functions require GeV
      // So make sure you pass them in GeV.
      double getWeight(const double x, const double y, const double Enu, const int targetZ); //in GeV
      //Initializer
      void read(enum modelSet m);

    private:
      double getWeightInternal(const double x, const double y, const double Enu, const int targetZ); //in GeV
  };

}

#endif

