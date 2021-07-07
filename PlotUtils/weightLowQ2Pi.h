#ifndef weightLowQ2Pi_h
#define weightLowQ2Pi_h

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "TGraph.h"
#include "TFile.h"


//For Compatibility with ROOT compiler uncomment the following:
//#define ROOT

namespace PlotUtils{

  //Class for getting the low Q2 suppression weights
  class weightLowQ2Pi
  {
    public:
      //Constructor: Read weights and weight shifts TGraphs from root file
      weightLowQ2Pi();


      // Provide Q2 in ()GeV/c)^2
          double getWeight(const double Q2, std::string channel, int shift = 0);
      double getMinosWeight(const double Q2 /*, int variation*/); //TODO MINOS uncertainties
                                                                  //TODO don't hardcode minos params

      //double Aminos, Aminos_err, Q0minos, Q0minos_err;

    private:
      TFile* weights_file;
      TGraph *h_cvweight_JOINT;
      TGraph *h_cvweight_shift_JOINT;
      TGraph *h_cvweight_NU1PI;
      TGraph *h_cvweight_shift_NU1PI;
      TGraph *h_cvweight_NUNPI;
      TGraph *h_cvweight_shift_NUNPI;
      TGraph *h_cvweight_NUPI0;
      TGraph *h_cvweight_shift_NUPI0;
      TGraph *h_cvweight_NUBARPI0;
      TGraph *h_cvweight_shift_NUBARPI0;

      void read(const std::string f);
  };

  // Static instance of lowq2pi weighter
  PlotUtils::weightLowQ2Pi& weight_lowq2pi();

}


#endif
