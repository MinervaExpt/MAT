#include "weight_2p2h.h"
#include "TSystem.h"

using namespace PlotUtils;

void weight_2p2h::read(const std::string f)
//Read in the params doubles from a file
//argument: valid filename
{
  double temp;
  int counter = 0;
  double params[6];
  
  TString filename = f + ".txt";  // expands environmentals
 
  std::cout << "weight_2p2h: 2p2h input file " << filename << std::endl;
  std::ifstream fh(filename.Data());
  if ( fh.is_open()){
    std::cout << "weight_2p2h: Reading from file " << f << std::endl;
    
    //Gather the params then assign afterward
    while (fh >> temp && counter < 6){
      params[counter] = temp;
      counter++;
    }
    assert(counter > 1);
    
    //Set members according to read
    norm = params[0];
    meanq0 = params[1];
    meanq3 = params[2];
    sigmaq0 = params[3];
    sigmaq3 = params[4];
    corr = params[5];
  }
  else{
    //File could not be read
    std::cout << "weight_2p2h: File could not be read" << std::endl;
    
  }
}


double weight_2p2h::Gaussian2D(const double q0, const double q3)
//FIXME: Compute the 2D Gaussian parameters ???
{
  //Unpack stats values from input array and combine w class members
  
  double z = (q0 - meanq0)*(q0 - meanq0) /sigmaq0/sigmaq0
  + (q3 - meanq3)*(q3 - meanq3) / sigmaq3/sigmaq3
  - 2*corr*(q0-meanq0)*(q3-meanq3)/ (sigmaq0 * sigmaq3);
  
  double ret = norm*exp( -0.5 * z / (1 - corr*corr) );
  
  //Need to add 1 to the results
  return ret;
}

// Static instances of 2p2h weighters
PlotUtils::weight_2p2h& PlotUtils::weight_2p2h_cv() {
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string directory_data = std::string(mparalocation)+"/data/Reweight/";
  static PlotUtils::weight_2p2h* _weight_2p2h_cv = 
    new PlotUtils::weight_2p2h(directory_data+"/fit-mec-2d-noScaleDown-penalty00300-best-fit");
  return *_weight_2p2h_cv;
}


PlotUtils::weight_2p2h& PlotUtils::weight_2p2h_nn() {
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string directory_data = std::string(mparalocation)+"/data/Reweight/";
  static PlotUtils::weight_2p2h* _weight_2p2h_nn = 
    new PlotUtils::weight_2p2h(directory_data+"/fit-mec-2d-nn-only-noScaleDown-penalty00300-best-fit");
  return *_weight_2p2h_nn;
}


PlotUtils::weight_2p2h& PlotUtils::weight_2p2h_np() {
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string directory_data = std::string(mparalocation)+"/data/Reweight/";
  static PlotUtils::weight_2p2h* _weight_2p2h_np = 
    new PlotUtils::weight_2p2h(directory_data+"/fit-mec-2d-np-only-noScaleDown-penalty02000-best-fit");
  return *_weight_2p2h_np;
}


PlotUtils::weight_2p2h& PlotUtils::weight_2p2h_qe() {
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string directory_data = std::string(mparalocation)+"/data/Reweight/";
  static PlotUtils::weight_2p2h* _weight_2p2h_qe = 
    new PlotUtils::weight_2p2h(directory_data+"/fit-qe-gaussian-noScaleDown-penalty02000-best-fit");
  return *_weight_2p2h_qe;
}
