#ifndef weight_minervaq2qe_h
#define weight_minervaq2qe_h

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
  This is a class for taking a q2qe value and reweighting values based on the 2D ME q2qe data/mc ratio.
 */

namespace PlotUtils{

  class weight_minervaq2qe
  {
    public:
       
      //Constructor: Read in params from a filename
    weight_minervaq2qe(const TString  f) { read(f); }    //Read in params from files
      
    TFile* fQ2QERatio;
    TH1D *hQ2QERatio;
    
    
    // MINERvA holds kinematics in MeV, but all these functions require GeV
    // So make sure you pass them in GeV.
    double getWeight(const double q2qe); //in GeV2
    //Initializer
    void read(TString filename);
    
  private:
    double getWeightInternal(const double q2qe); //in GeV2
  };
  
}

#endif

