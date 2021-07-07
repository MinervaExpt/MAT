#ifndef MichelSYSTEMATICS_H
#define MichelSYSTEMATICS_H

#include "TreeWrapper.h"
#include "NSFDefaults.h"

namespace PlotUtils{

  template<class T>
  class MichelEfficiencyUniverse: public T
    {
    public:
      MichelEfficiencyUniverse(typename T::config_t chw, double nsigma, double fracUncertainty = NSFDefaults::MichelTagEfficiency_Err );
      
      virtual double GetMichelEfficiencyWeight() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return true; }/*override*/
    private:
      double m_fracUncertainty;
    };
  
}

// Explicit Instantiation, needed for loading into python
//! Make sure you put this into Ana/PlotUtils/dict/PlotUtilsDict.h
//template class

// For template classes, the header needs to know about the definitions
#include "MichelSystematics.cxx"

#endif // AngleSYSTEMATICS_H

