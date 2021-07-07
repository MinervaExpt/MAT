#ifndef MUONSYSTEMATICS_H
#define MUONSYSTEMATICS_H

#include "TreeWrapper.h"
#include "NSFDefaults.h"
#include "FluxSystematics.cxx"

// Helper functions declared/defined in the .cxx
// GetMinervaMuonSystematicsMap( typename T::config_t chain );
// GetMinosMuonSystematicsMap( typename T::config_t chain );))

namespace PlotUtils{
  //=================================================================================
  // Minerva muon-momentum-shifted universe 
  //=================================================================================
  template<class T>
  class MuonUniverseMinerva: public T
  {
    public:
      MuonUniverseMinerva(typename T::config_t chw, double nsigma );

      double GetMuonMomentumShiftMinerva() const;
      virtual double GetPmuMinerva() const /*override*/;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
  };
  //=================================================================================
  // Minos muon-momentum-shifted universe
  //=================================================================================
  template<class T>
  class MuonUniverseMinos: public T
  {
    public:
      MuonUniverseMinos(typename T::config_t chw, double nsigma );

      double GetMuonMomentumShiftMinos() const;
      virtual double GetFluxAndCVWeight( double Enu = -99. /*GeV*/, int nu_pdg = -99) const /* override */; 
      virtual double GetPmuMinos() const /*override*/;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "MuonSystematics.cxx"

#endif // MUONSYSTEMATICS_H

