#ifndef TARGETMASSSYSTEMATICS_H
#define TARGETMASSSYSTEMATICS_H

#include "TreeWrapper.h"

// Helper functions declared/defined in the .cxx
// GetTargetMassSystematicsMap(typename T::config_t chain );

namespace PlotUtils{
  template<class T>
  class TargetMassUniverse: public T
  {
    public:
      TargetMassUniverse(typename T::config_t chw, double nsigma );

      virtual double GetTargetMassWeight() const /*override*/;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;
  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "TargetMassSystematics.cxx"

#endif // TARGETMASSSYSTEMATICS_H

