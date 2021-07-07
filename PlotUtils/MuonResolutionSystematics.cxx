#ifndef MuonResolutionSYSTEMATICS_CXX
#define MuonResolutionSYSTEMATICS_CXX

#include "MuonResolutionSystematics.h"
#include "TVector3.h" // Using this to get theta from a 3d vector
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetMuonResolutionSystematics(typename T::config_t chain,  double fracUncertainty=NSFDefaults::muonResolution_Err) {
    std::vector<T*> ret;
    std::cout << "Muon Resolution Systematics created with fractional uncertainty " << fracUncertainty << std::endl;
    ret.push_back(new PlotUtils::MuonResolutionUniverse<T>(chain, -1., fracUncertainty));
    ret.push_back(new PlotUtils::MuonResolutionUniverse<T>(chain, 1., fracUncertainty));
    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMuonResolutionSystematicsMap(typename T::config_t chain, double fracUncertainty=NSFDefaults::muonResolution_Err) {

    std::cout << "Muon Resolution Systematics created with fractional uncertainty " << fracUncertainty << std::endl;
    std::map< std::string, std::vector<T*> > ret;
    
    ret["MuonResolution"].push_back(new PlotUtils::MuonResolutionUniverse<T>(chain, -1., fracUncertainty));
    ret["MuonResolution"].push_back(new PlotUtils::MuonResolutionUniverse<T>(chain, 1., fracUncertainty));

    return ret;
  }

  template <class T>
  std::vector<T*> GetMuonAngleResolutionSystematics(typename T::config_t chain, double fracUncertainty=NSFDefaults::muon_angle_frac_res) {
    std::vector<T*> ret;
    ret.push_back(new PlotUtils::MuonAngleXResolutionUniverse<T>(chain, -1.,fracUncertainty));
    ret.push_back(new PlotUtils::MuonAngleXResolutionUniverse<T>(chain, 1., fracUncertainty));
    ret.push_back(new PlotUtils::MuonAngleYResolutionUniverse<T>(chain, -1.,fracUncertainty));
    ret.push_back(new PlotUtils::MuonAngleYResolutionUniverse<T>(chain, 1., fracUncertainty));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMuonAngleResolutionSystematicsMap(typename T::config_t chain,  double fracUncertainty=NSFDefaults::muon_angle_frac_res) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["MuonAngleXResolution"].push_back(new PlotUtils::MuonAngleXResolutionUniverse<T>(chain, -1., fracUncertainty));
    ret["MuonAngleXResolution"].push_back(new PlotUtils::MuonAngleXResolutionUniverse<T>(chain, 1.,  fracUncertainty));
    ret["MuonAngleYResolution"].push_back(new PlotUtils::MuonAngleYResolutionUniverse<T>(chain, -1., fracUncertainty));
    ret["MuonAngleYResolution"].push_back(new PlotUtils::MuonAngleYResolutionUniverse<T>(chain, 1.,  fracUncertainty));

    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::MuonResolutionUniverse<T>::MuonResolutionUniverse(
    typename T::config_t chw, double nsigma, double fracUncertainty) : T(chw, nsigma), m_fracUncertainty(fracUncertainty) {}


template<typename T>
double PlotUtils::MuonResolutionUniverse<T>::GetMuonResolutionMomentumShift() const {
  double modP;
  // Get shift due to muon resolution
  double recoP = T::GetPmu();
  std::vector<double> truePVec = T::GetVecDouble("mc_primFSLepton");

  double trueP = std::sqrt(truePVec[0]*truePVec[0]+ truePVec[1]*truePVec[1]+truePVec[2]*truePVec[2]);

  modP = (trueP-recoP)*(T::m_nsigma*m_fracUncertainty);

  double shift_val = modP;
  
  return shift_val;
}

template<typename T>
double PlotUtils::MuonResolutionUniverse<T>::GetPmu() const {

  double shift_val = GetMuonResolutionMomentumShift();

  return shift_val+T::GetPmu();
}

template<typename T>
std::string PlotUtils::MuonResolutionUniverse<T>::ShortName() const { return "Muon_Energy_Resolution"; }

template<typename T>
std::string PlotUtils::MuonResolutionUniverse<T>::LatexName() const { return "Muon_Energy_Resolution"; }

template<typename T>
PlotUtils::MuonAngleXResolutionUniverse<T>::MuonAngleXResolutionUniverse(
    typename T::config_t chw, double nsigma, double fracUncertainty) : T(chw, nsigma), m_fracUncertainty(fracUncertainty) {}

template<typename T>
double PlotUtils::MuonAngleXResolutionUniverse<T>::GetTrueThetaXmu() const {
  std::vector<double> truePVec = T::GetVecDouble("mc_primFSLepton");
  TVector3 p3mu( truePVec[0], truePVec[1], truePVec[2] );
  p3mu.RotateX(MinervaUnits::numi_beam_angle_rad);

  double denom2 = p3mu.X()*p3mu.X() + p3mu.Z()*p3mu.Z();
  if( denom2 == 0 ) return -999; //when we don't have any momentum

  double trueThetaX = p3mu.X() > 0 ?  std::acos( p3mu.Z()/ sqrt(denom2) ) :
                                     -std::acos( p3mu.Z()/ sqrt(denom2) );
  return trueThetaX;
}

template<typename T>
double PlotUtils::MuonAngleXResolutionUniverse<T>::GetThetaXmuShift() const {

  double recoThetaX = T::GetThetaXmu();
  double trueThetaX = GetTrueThetaXmu();
  double modThetaX = (trueThetaX - recoThetaX)*(T::m_nsigma*m_fracUncertainty); 

  return modThetaX;
}

template<typename T>
double PlotUtils::MuonAngleXResolutionUniverse<T>::GetThetaXmu() const {
  double shift_val = GetThetaXmuShift();

  return shift_val+T::GetThetaXmu();
}

template<typename T>
std::string PlotUtils::MuonAngleXResolutionUniverse<T>::ShortName() const { return "MuonAngleXResolution"; }

template<typename T>
std::string PlotUtils::MuonAngleXResolutionUniverse<T>::LatexName() const { return "Muon Track Angle X Resolution (rad.)"; }

template<typename T>
PlotUtils::MuonAngleYResolutionUniverse<T>::MuonAngleYResolutionUniverse(
    typename T::config_t chw, double nsigma, double fracUncertainty) : T(chw, nsigma), m_fracUncertainty(fracUncertainty) {}

template<typename T>
double PlotUtils::MuonAngleYResolutionUniverse<T>::GetTrueThetaYmu() const {
  std::vector<double> truePVec = T::GetVecDouble("mc_primFSLepton");
  TVector3 p3mu( truePVec[0], truePVec[1], truePVec[2] );
  p3mu.RotateX(MinervaUnits::numi_beam_angle_rad);

  double denom2 = p3mu.Y()*p3mu.Y() + p3mu.Z()*p3mu.Z();
  if( denom2 == 0 ) return -999; //when we don't have any momentum

  double trueThetaY = p3mu.Y() > 0 ?  std::acos( p3mu.Z()/ sqrt(denom2) ) :
                                     -std::acos( p3mu.Z()/ sqrt(denom2) );
  return trueThetaY;
}

template<typename T>
double PlotUtils::MuonAngleYResolutionUniverse<T>::GetThetaYmuShift() const {

  double recoThetaY = T::GetThetaYmu();
  double trueThetaY = GetTrueThetaYmu();
  double modThetaY = (trueThetaY - recoThetaY)*(T::m_nsigma*m_fracUncertainty); 

  return modThetaY;
}

template<typename T>
double PlotUtils::MuonAngleYResolutionUniverse<T>::GetThetaYmu() const {
  double shift_val = GetThetaYmuShift();

  return shift_val+T::GetThetaYmu();
}

template<typename T>
std::string PlotUtils::MuonAngleYResolutionUniverse<T>::ShortName() const { return "MuonAngleYResolution"; }

template<typename T>
std::string PlotUtils::MuonAngleYResolutionUniverse<T>::LatexName() const { return "Muon Track Angle X Resolution (rad.)"; }


#endif // MuonResolutionSYSTEMATICS_CXX
