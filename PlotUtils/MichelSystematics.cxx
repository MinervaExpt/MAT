#ifndef MichelSystematics_CXX
#define MichelSystematics_CXX

#include "MichelSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetMichelEfficiencySystematics(typename T::config_t chain, double fracUncertainty = NSFDefaults::MichelTagEfficiency_Err ) {
    std::vector<T*> ret;
    std::cout << "Michel Tag Efficiency Uncertainty " <<  fracUncertainty << "\%" << std::endl;
    ret.push_back(new PlotUtils::MichelEfficiencyUniverse<T>(chain, -1., fracUncertainty));
    ret.push_back(new PlotUtils::MichelEfficiencyUniverse<T>(chain, 1., fracUncertainty));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMichelEfficiencySystematicsMap(typename T::config_t chain, double fracUncertainty = NSFDefaults::MichelTagEfficiency_Err ) {
    std::map< std::string, std::vector<T*> > ret;
    std::cout << "Michel Tag Efficiency Uncertainty " <<  fracUncertainty << "\%" << std::endl;
    
    ret["MichelEfficiency"].push_back(new PlotUtils::MichelEfficiencyUniverse<T>(chain, -1., fracUncertainty));
    ret["MichelEfficiency"].push_back(new PlotUtils::MichelEfficiencyUniverse<T>(chain, 1., fracUncertainty));

    return ret;
  }

}

//!  Michel Tag Efficiency
// Class Definitions
// Constructor
template<typename T>
PlotUtils::MichelEfficiencyUniverse<T>::MichelEfficiencyUniverse(
    typename T::config_t chw, double nsigma, double fracUncertainty) : T(chw, nsigma), m_fracUncertainty(fracUncertainty) {}

//Note, there is no EM shifts here because that was already factored into the efficiency calculation
//There was a 1.5% uncertainty on the michel energy
template<typename T>
double PlotUtils::MichelEfficiencyUniverse<T>::GetMichelEfficiencyWeight() const {
  if(T::IsTruth()) return 1; //This is a detector effect only
  double shift_val = 1 + T::m_nsigma * m_fracUncertainty;//probably the easiest reweight ever
  return shift_val;
}

template<typename T>
std::string PlotUtils::MichelEfficiencyUniverse<T>::ShortName() const { return "MichelEfficiency"; }

template<typename T>
std::string PlotUtils::MichelEfficiencyUniverse<T>::LatexName() const { return "Michel Tag Efficiency"; }

#endif // MichelSystematics_CXX
