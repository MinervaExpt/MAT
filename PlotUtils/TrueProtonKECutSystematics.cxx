#ifndef TrueProtonKECutSystematics_CXX
#define TrueProtonKECutSystematics_CXX

#include "TrueProtonKECutSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetTrueProtonKECutSystematics(typename T::config_t chain,
                                     double Uncertainty =NSFDefaults::TrueProtonKECutVariation ) {
    std::vector<T*> ret;
    std::cout << "True Proton KE Cut Systematics created with CV cut at" << T::GetTrueProtonKECutCentral() << " and uncertainty "  <<  Uncertainty << " MeV" << std::endl;
    ret.push_back(new PlotUtils::TrueProtonKECutUniverse<T>(chain, -1., Uncertainty));
    ret.push_back(new PlotUtils::TrueProtonKECutUniverse<T>(chain, 1., Uncertainty));
     return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetTrueProtonKECutSystematicsMap(
      typename T::config_t chain, double Uncertainty=NSFDefaults::TrueProtonKECutVariation) {
    std::map< std::string, std::vector<T*> > ret;
    std::cout << "True Proton KE Cut Systematics created with CV cut at" << T::GetTrueProtonKECutCentral() << " and uncertainty "  << Uncertainty  <<  " MeV" <<   std::endl;
    
    ret["TrueProtonKECut"].push_back(new PlotUtils::TrueProtonKECutUniverse<T>(chain, -1., Uncertainty));
    ret["TrueProtonKECut"].push_back(new PlotUtils::TrueProtonKECutUniverse<T>(chain, 1., Uncertainty));
    
    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::TrueProtonKECutUniverse<T>::TrueProtonKECutUniverse(
    typename T::config_t chw, double nsigma, double Uncertainty) : T(chw, nsigma), m_Uncertainty(Uncertainty) {}

template<typename T>
double PlotUtils::TrueProtonKECutUniverse<T>::GetTrueProtonKECut() const {

  double shift_val = T::m_nsigma * m_Uncertainty;
  //std::cout << "shift proton" << shift_val << std::endl;
  return shift_val+T::GetTrueProtonKECut();
  
}

template<typename T>
std::string PlotUtils::TrueProtonKECutUniverse<T>::ShortName() const { return "TrueProtonKECut"; }


template<typename T>
std::string PlotUtils::TrueProtonKECutUniverse<T>::LatexName() const { return "True Proton KE Cut, MeV"; }



#endif // TrueProtonKECutSystematics_CXX
