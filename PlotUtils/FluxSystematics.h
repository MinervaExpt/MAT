#ifndef FLUXSYSTEMATICS_H
#define FLUXSYSTEMATICS_H

#include "TreeWrapper.h"
#include "FluxReweighter.h"

// Helper functions declared/defined in the .cxx
// GetFluxSystematicsMap(typename T::config_t chain, unsigned int n_universes);

namespace PlotUtils{

  template<class T>
  class FluxUniverse : public T
  {
    public:
      FluxUniverse(typename T::config_t chw, double nsigma, int universe_number = -1);

    FluxUniverse(){ //std::cout << " make a Flux Universe " << std::endl;
    }
      // I can't get things to compile with "override"
      virtual double GetFluxAndCVWeight( double Enu = -99. /*GeV*/, int nu_pdg = -99) const /* override */; 

      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      void SetUniverseNumber(int i);
      int m_universe_number;
      bool m_use_nuE_constraint_error;
  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "FluxSystematics.cxx"

#endif // FLUXSYSTEMATICS_H
