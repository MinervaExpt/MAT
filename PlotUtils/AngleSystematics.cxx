#ifndef AngleSystematics_CXX
#define AngleSystematics_CXX

#include "AngleSystematics.h"
#include "TRandom3.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetAngleSystematics(typename T::config_t chain,
                                     double XUncertainty=NSFDefaults::beamThetaX_Err, double YUncertainty=NSFDefaults::beamThetaY_Err) {
    std::vector<T*> ret;
std::cout << "Angle Systematics created with thetaX err " <<  XUncertainty << " and thetaY err " <<  YUncertainty <<  std::endl;
    ret.push_back(new PlotUtils::BeamAngleXUniverse<T>(chain, -1., XUncertainty));
    ret.push_back(new PlotUtils::BeamAngleXUniverse<T>(chain, 1., XUncertainty));
    ret.push_back(new PlotUtils::BeamAngleYUniverse<T>(chain, -1., YUncertainty));
    ret.push_back(new PlotUtils::BeamAngleYUniverse<T>(chain, 1., YUncertainty));


    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetAngleSystematicsMap(
      typename T::config_t chain, double XUncertainty=NSFDefaults::beamThetaX_Err, double YUncertainty=NSFDefaults::beamThetaY_Err) {
    std::map< std::string, std::vector<T*> > ret;
    std::cout << "Angle Systematics created with thetaX err " <<  XUncertainty << " and thetaY err " <<  YUncertainty <<  std::endl;
    
    ret["BeamAngleX"].push_back(new PlotUtils::BeamAngleXUniverse<T>(chain, -1., XUncertainty));
    ret["BeamAngleX"].push_back(new PlotUtils::BeamAngleXUniverse<T>(chain, 1., XUncertainty));
    ret["BeamAngleY"].push_back(new PlotUtils::BeamAngleYUniverse<T>(chain, -1., YUncertainty));
    ret["BeamAngleY"].push_back(new PlotUtils::BeamAngleYUniverse<T>(chain, 1., YUncertainty));

    return ret;
  }

}

//!  BeamAngleOffset X
// Class Definitions
// Constructor
template<typename T>
PlotUtils::BeamAngleXUniverse<T>::BeamAngleXUniverse(
    typename T::config_t chw, double nsigma, double Uncertainty) : T(chw, nsigma), m_Uncertainty(Uncertainty) {}

template<typename T>
double PlotUtils::BeamAngleXUniverse<T>::GetBeamAngleOffsetX() const {

  double shift_val = T::m_nsigma * m_Uncertainty;

  return shift_val+T::GetBeamAngleOffsetX();
}

template<typename T>
std::string PlotUtils::BeamAngleXUniverse<T>::ShortName() const { return "BeamAngleX"; }

template<typename T>
std::string PlotUtils::BeamAngleXUniverse<T>::LatexName() const { return "Beam Angle X (rad.)"; }


//!  AngleY

template<typename T>
PlotUtils::BeamAngleYUniverse<T>::BeamAngleYUniverse(
    typename T::config_t chw, double nsigma, double Uncertainty) : T(chw, nsigma), m_Uncertainty(Uncertainty) {}

template<typename T>
double PlotUtils::BeamAngleYUniverse<T>::GetBeamAngleOffsetY() const {

  double shift_val = T::m_nsigma * m_Uncertainty;

  return shift_val+T::GetBeamAngleOffsetY();
}

template<typename T>
std::string PlotUtils::BeamAngleYUniverse<T>::ShortName() const { return "BeamAngleY"; }

template<typename T>
std::string PlotUtils::BeamAngleYUniverse<T>::LatexName() const { return "Beam Angle Y (rad.)"; }

#endif // AngleSystematics_CXX
