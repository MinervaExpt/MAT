#ifndef weightNuclearScreening_h
#define weightNuclearScreening_h

#include <fstream>  //ifstream
#include <iostream> //cout
//#include <cstdlib>  //cout, endl

#include <string>
#include "math.h"
#include "assert.h"
#include "TFile.h"
#include "TH2D.h"

//For Compatibility with ROOT compiler uncomment the following:
//#define ROOT


//Class for getting Gaussian paramaters inputs from a given file

namespace PlotUtils{
  
  class weightNuclearScreening{
  public:

    //Constructor: Read in params from a filename
    weightNuclearScreening(const std::string f){read(f);}    
    void read(const std::string f);
    TFile* fNucScrFile;
    TH2D* fNucScrRatio;
    double GetWeight(double nu, double Q2);
  };
}

#endif

