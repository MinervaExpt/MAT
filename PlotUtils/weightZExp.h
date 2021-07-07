#ifndef weightZExp_h
#define weightZExp_h

#include <fstream>  //ifstream
#include <iostream> //cout

#include <TString.h>
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
  This is a class for taking a Q2 histogram(s) and applying a variety of weights to reflex the transformation of 
CCQE with dipole parameterization of the axial FF and the z-expansion form via 
Bhubanjyoti Bhattacharya, Richard J. Hill, and Gil Paz Phys. Rev. D 84, 073006 – Published 13 October 2011
Aaron S. Meyer, Minerba Betancourt, Richard Gran, and Richard J. Hill Phys. Rev. D 93, 113015 – Published 23 June 2016
 */

namespace PlotUtils{

  class weightZExp
  {
    public:
      //Constructor: Read in params from a filename
      weightZExp(const TString  f) { read(f); }    //Read in params from file
      
      TString  filename;
      TFile* fZExpRatioCov;
      PlotUtils::MnvH1D *hZExpRatioCov;
      PlotUtils::MnvVertErrorBand *hZExpRatioBand;

      // MINERvA holds kinematics in MeV, but all these functions require GeV
      // So make sure you pass them in GeV.
      double getWeight(const double mc_q2); //in GeV
      double getWeight(const double mc_q2, const int uni_n); //in GeV

      //Initializer
      void read(const TString  f);

    private:
      double getWeightInternal(const double mc_Q2, const int uni_n);
      double getWeightInternal(const double mc_Q2);
  };

}

#endif

