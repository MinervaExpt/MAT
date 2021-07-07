#ifndef TARGETMASSSYSTEMATICS_CXX
#define TARGETMASSSYSTEMATICS_CXX

#include "NSFDefaults.h"
#include "TargetMassSystematics.h"
#include "TargetUtils.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetTargetMassSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;
    
    ret.push_back(new PlotUtils::TargetMassUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassUniverse<T>(chain, 1.));

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetTargetMassSystematicsMap(
      typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["TargetMass"].push_back(new PlotUtils::TargetMassUniverse<T>(chain, -1.));
    ret["TargetMass"].push_back(new PlotUtils::TargetMassUniverse<T>(chain, 1.));

    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::TargetMassUniverse<T>::TargetMassUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}


template<typename T>
double PlotUtils::TargetMassUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus
  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // if interaction is not in Nuclear Target regio, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // if interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if( PlotUtils::TargetUtils::Get().InWaterTargetMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                         T::GetVecElem("mc_vtx",2), targetZ ) ){
    wgt_shift = NSFDefaults::h2o_err;
  }
  else if(targetZ==6){
    wgt_shift = NSFDefaults::c_err; 
  }
  else if(targetZ==26){
    wgt_shift = NSFDefaults::fe_err; 
  }
  else if(targetZ==82){
    wgt_shift = NSFDefaults::pb_err; 
  }
  // assume CH if not H2O, C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

//This systematic is special: its weight function is always 1 in the CV
//as of April 2021.  So, it doesn't even have a Reweighter.  Just
//including it in the list of systematics and using PlotUtils::Model applies it.
template <typename T>
double PlotUtils::TargetMassUniverse<T>::GetWeightRatioToCV() const {
  double cv_wgt = 1; //TODO: Watch very carefully for this to change

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus
  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // if interaction is not in Nuclear Target regio, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // if interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if(targetZ==6){
    wgt_shift = NSFDefaults::c_err;
  }
  else if(targetZ==26){
    wgt_shift = NSFDefaults::fe_err;
  }
  else if(targetZ==82){
    wgt_shift = NSFDefaults::pb_err;
  }
  // assume CH if not C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  1 + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassUniverse<T>::ShortName() const { return "Target_Mass"; }


template<typename T>
std::string PlotUtils::TargetMassUniverse<T>::LatexName() const { return "Target Mass"; }


#endif // TARGETMASSSYSTEMATICS_CXX
