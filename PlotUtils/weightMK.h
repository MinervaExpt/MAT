#ifndef weightMK_h
#define weightMK_h

#include <iostream>
#include <fstream>  //ifstream
#include <stdexcept>
#include <cassert>

#include <TString.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "math.h"
#include "assert.h"

#include "TGraph.h"
#include "TFile.h"
#include <Rtypes.h>
#include <vector>

namespace PlotUtils {

  class TreeWrapper;

  //Class for getting the ratios of MK model/GENIE
  //The Weight is based in Single pion production in neutrino-nucleon interactions
  //M. Kabirnezhad
  //Phys. Rev. D 97, 013002 – Published 23 January 2018
  //https://doi.org/10.1103/PhysRevD.97.013002
  class weightMK
  {
    public:
      //Constructor: Read in params from a filename
      weightMK(const TString  f) { read(f) ;}
      
      TString  filename;
      TFile* fMKratio;

      TH2F *hMKratio_CC1ppip;
      TH2F *hMKratio_CC1pi0;
      TH2F *hMKratio_CC1npip;

      // Make sure you pass them in GeV.
      // and also mc_W and mc_Q2 are GENIE true values

      double getWeight(const double mc_W/*should be true GENIE mc_w */,
		       const double mc_Q2/*should be true GENIE mc_Q2 */,
		       int variation); // in GeV

      double getWeight(const double mc_W/*should be true GENIE mc_w */,
		       const double mc_Q2/*should be true GENIE mc_Q2 */, 
		       const int mc_current, 
		       const int mc_intType,
		       const int mc_er_nPart, 
		       const int mc_targetA, 
		       const int *mc_er_ID, 
		       const int *mc_er_status); // in GeV

      double getWeight(const double mc_W/*should be true GENIE mc_w */,
		       const double mc_Q2/*should be true GENIE mc_Q2 */, 
		       const int mc_current, 
		       const int mc_intType,
		       const int mc_er_nPart, 
		       const int mc_targetA, 
		       const std::vector<int>& mc_er_ID, 
		       const std::vector<int>& mc_er_status); // in GeV

      double getWeight(PlotUtils::TreeWrapper* chw,
		       Long64_t entry);



      void read(const TString f);

    private:
      double getWeightInternal(const double mc_W, const double mc_Q2, int variation);

  };

  // Static instance of MK weighter
  PlotUtils::weightMK& weight_mk();

}

#endif
