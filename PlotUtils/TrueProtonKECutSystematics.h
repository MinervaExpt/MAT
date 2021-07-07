#ifndef TrueProtonKECutSYSTEMATICS_H
#define TrueProtonKECutSYSTEMATICS_H

// class to vary the TrueProtonKECut

#include "TreeWrapper.h"
#include "NSFDefaults.h"

namespace PlotUtils{
  template<class T>
  class TrueProtonKECutUniverse: public T
    {
    public:
      TrueProtonKECutUniverse(typename T::config_t chw, double nsigma , double Uncertainty=NSFDefaults::TrueProtonKECutVariation);
      virtual double GetTrueProtonKECut() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_Uncertainty;
    };
}


// For template classes, the header needs to know about the definitions
#include "TrueProtonKECutSystematics.cxx"

#endif // TrueProtonKECutSYSTEMATICS_H

