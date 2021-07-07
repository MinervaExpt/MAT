#ifndef RecoProtonKECutSYSTEMATICS_H
#define RecoProtonKECutSYSTEMATICS_H

// class to vary the RecoProtonKECut

#include "TreeWrapper.h"
#include "NSFDefaults.h"

namespace PlotUtils{
  template<class T>
  class RecoProtonKECutUniverse: public T
    {
    public:
      RecoProtonKECutUniverse(typename T::config_t chw, double nsigma , double Uncertainty=NSFDefaults::RecoProtonKECutVariation);
      virtual double GetRecoProtonKECut() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_Uncertainty;
    };
}


// For template classes, the header needs to know about the definitions
#include "RecoProtonKECutSystematics.cxx"

#endif // RecoProtonKECutSYSTEMATICS_H

