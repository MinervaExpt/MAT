#ifndef weight_2p2h_h
#define weight_2p2h_h

#include <fstream>  //ifstream
#include <iostream> //cout
//#include <cstdlib>  //cout, endl

#include <string>
#include "math.h"
#include "assert.h"

//For Compatibility with ROOT compiler uncomment the following:
//#define ROOT


//Class for getting Gaussian paramaters inputs from a given file

namespace PlotUtils{

  class weight_2p2h {
  public:
     //Constructor: Read in params from a filename
     weight_2p2h(const std::string f) { read(f); }    //Read in params from file

     std::string filename;
     double norm, meanq0, meanq3;
     double sigmaq0, sigmaq3, corr;

     //Initializer
     void read(const std::string f);

     // This is the 2D Gaussian weight function. The x[] array should be:
     // x[0] = true q3 in GeV
     // x[1] = true q0 in GeV
     //
     // The params[] array should be the best-fit parameters from one of
     // the fits. The order is as shown at the start of the function
     double Gaussian2D(const double q0, const double q3);
    
    double getWeight(const double q0, const double q3){return Gaussian2Dplusone(q0,q3);}

     // The 2d gaussian +1: the function used for the "no scale down" fits
     double Gaussian2Dplusone(const double q0, const double q3) {return 1. + Gaussian2D(q0,q3);}
  };

  // Static instances of 2p2h weighters
  weight_2p2h& weight_2p2h_cv();
  weight_2p2h& weight_2p2h_nn();
  weight_2p2h& weight_2p2h_np();
  weight_2p2h& weight_2p2h_qe();

}

#endif

