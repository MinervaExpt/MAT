#ifndef FLUXSYSTEMATICS_CXX
#define FLUXSYSTEMATICS_CXX

#include "PlotUtils/FluxSystematics.h"
#include <iostream>

// Helper Functions

namespace PlotUtils{
  // Get a flux reweighter for each playlist.
  // Technically, we only need one FRW for each "era", but eras could change.
  
  template <class T>
  std::vector<T*> GetFluxSystematics(typename T::config_t chain, 
                                     unsigned int n_universes){
    std::vector<T*> ret;
    double nsigma = 1.;
    for (unsigned int i = 0; i < n_universes; ++i)
      ret.push_back(new PlotUtils::FluxUniverse<T>(chain, nsigma, i));
    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetFluxSystematicsMap(typename T::config_t chain, 
                                                                 unsigned int n_universes) {
    std::map< std::string, std::vector<T*> > ret;
    const double nsigma = 1.;
    for (unsigned int i = 0; i < n_universes; ++i)
      ret["Flux"].push_back(new PlotUtils::FluxUniverse<T>(chain, nsigma, i));
    return ret;
  }

}


// Constuctor
template <class T>
PlotUtils::FluxUniverse<T>::FluxUniverse(typename T::config_t chw,
                                         double nsigma, int universe_number )
  : T(chw, nsigma), 
    m_universe_number(universe_number)
{//std::cout << " make a PlotUtils Flux Universe " << std::endl;
}


// Get Flux Weight -- pass plist string and nuE constraint on/off bool
template <class T>
double PlotUtils::FluxUniverse<T>::GetFluxAndCVWeight(double Enu,
                                                      int nu_pdg) const {
  if(Enu == -99.) Enu = T::GetDouble("mc_incomingE")*1e-3; // If the user hasn't specified enu, use default
  if(nu_pdg == -99) nu_pdg = T::GetInt("mc_incoming");     // If the user hasn't specific nu_pdg, use default
  if(m_universe_number < 0) 
    std::cout << "GetFluxAndCVWeight WARNING: universe number not "
                 "set.\nThis should have been done in the FluxUniverse "
                 "constructor or with FluxUniverse::SetUniverseNumber.\n";
  
  // FRW::GetFluxErrorWeight() returns a weight which is the ratio for a given Enu b/w the sys. universe wgt and the CV wgt
  // So, the correct usage is to multiply this weight by the CV weight, not to replace the CV weight with this one
  double sysWgt = PlotUtils::flux_reweighter(T::GetPlaylist(),nu_pdg,T::UseNuEConstraint(),T::GetNFluxUniverses()).GetFluxErrorWeight(Enu, nu_pdg, 
                                                                                   m_universe_number);
  return sysWgt*T::GetFluxAndCVWeight();
}

template <class T>
double PlotUtils::FluxUniverse<T>::GetWeightRatioToCV() const
{
  /*if(Enu == -99.)*/ const double Enu = T::GetDouble("mc_incomingE")*1e-3; // If the user hasn't specified enu, use default
  /*if(nu_pdg == -99)*/ const int nu_pdg = T::GetInt("mc_incoming");     // If the user hasn't specific nu_pdg, use default
  if(m_universe_number < 0)
    std::cout << "GetFluxAndCVWeight WARNING: universe number not "
                 "set.\nThis should have been done in the FluxUniverse "
                 "constructor or with FluxUniverse::SetUniverseNumber.\n";

  // FRW::GetFluxErrorWeight() returns a weight which is the ratio for a given Enu b/w the sys. universe wgt and the CV wgt
  // So, the correct usage is to multiply this weight by the CV weight, not to replace the CV weight with this one
  return PlotUtils::flux_reweighter(T::GetPlaylist(),nu_pdg,T::UseNuEConstraint(),T::GetNFluxUniverses()).GetFluxErrorWeight(Enu, nu_pdg,
                                                                                   m_universe_number);
}


template <class T>
void PlotUtils::FluxUniverse<T>::SetUniverseNumber(int i) {
  m_universe_number = i;
}


template <class T>
std::string PlotUtils::FluxUniverse<T>::ShortName() const {
  return "Flux";
}


template <class T>
std::string PlotUtils::FluxUniverse<T>::LatexName() const {
  return "Flux";
}

#endif // FLUXSYSTEMATICS_CXX
