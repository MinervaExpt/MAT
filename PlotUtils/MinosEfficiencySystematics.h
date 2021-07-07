#ifndef MINOSEFFICIENCYSYSTEMATICS_H
#define MINOSEFFICIENCYSYSTEMATICS_H

#include "TreeWrapper.h"

// Helper functions declared/defined in the .cxx
// GetMinosEfficiencySystematicsMap(typename T::config_t chain );

namespace PlotUtils{
  template<class T>
  class MinosEfficiencyUniverse: public T
  {
    public:
      MinosEfficiencyUniverse(typename T::config_t chw, double nsigma );

      virtual double GetMinosEfficiencyWeight() const /*override*/;
      virtual double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;
  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "MinosEfficiencySystematics.cxx"

#endif // MINOSEFFICIENCYSYSTEMATICS_H

