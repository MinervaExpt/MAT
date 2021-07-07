#ifndef RecoProtonKECutSystematics_CXX
#define RecoProtonKECutSystematics_CXX

#include "RecoProtonKECutSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetRecoProtonKECutSystematics(typename T::config_t chain,
                                     double Uncertainty =NSFDefaults::RecoProtonKECutVariation ) {
    std::vector<T*> ret;
    std::cout << "Reco Proton KE Cut Systematics created with CV cut at" << T::GetRecoProtonKECutCentral() << " and uncertainty "  <<  Uncertainty << " MeV" << std::endl;
    ret.push_back(new PlotUtils::RecoProtonKECutUniverse<T>(chain, -1., Uncertainty));
    ret.push_back(new PlotUtils::RecoProtonKECutUniverse<T>(chain, 1., Uncertainty));
     return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetRecoProtonKECutSystematicsMap(
      typename T::config_t chain, double Uncertainty=NSFDefaults::RecoProtonKECutVariation) {
    std::map< std::string, std::vector<T*> > ret;
    std::cout << "Reco Proton KE Cut Systematics created with CV cut at" << T::GetRecoProtonKECutCentral() << " and uncertainty "  << Uncertainty  <<  " MeV" <<   std::endl;
    
    ret["RecoProtonKECut"].push_back(new PlotUtils::RecoProtonKECutUniverse<T>(chain, -1., Uncertainty));
    ret["RecoProtonKECut"].push_back(new PlotUtils::RecoProtonKECutUniverse<T>(chain, 1., Uncertainty));
    
    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::RecoProtonKECutUniverse<T>::RecoProtonKECutUniverse(
    typename T::config_t chw, double nsigma, double Uncertainty) : T(chw, nsigma), m_Uncertainty(Uncertainty) {}

template<typename T>
double PlotUtils::RecoProtonKECutUniverse<T>::GetRecoProtonKECut() const {

  double shift_val = T::m_nsigma * m_Uncertainty;
  //std::cout << "shift proton" << shift_val << std::endl;
  return shift_val+T::GetRecoProtonKECut();
  
}

template<typename T>
std::string PlotUtils::RecoProtonKECutUniverse<T>::ShortName() const { return "RecoProtonKECut"; }


template<typename T>
std::string PlotUtils::RecoProtonKECutUniverse<T>::LatexName() const { return "Reco Proton KE Cut, MeV"; }



#endif // RecoProtonKECutSystematics_CXX
