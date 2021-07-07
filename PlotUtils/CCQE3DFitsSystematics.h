#ifndef SYSTEMATICS_3DCQEFITS_H
#define SYSTEMATICS_3DCQEFITS_H

#include "TreeWrapper.h"
#include <TH2D.h>

namespace PlotUtils{
  template<class T>
  class CCQE3DFitsUniverse: public T
  {
    public:
      CCQE3DFitsUniverse(typename T::config_t chw, double nsigma, int variation );

      virtual double GetCCQE3DFitsWeight() const /*override*/;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      int m_variation;

      std::vector<TH2D*> m_weightHists; 

  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "CCQE3DFitsSystematics.cxx"

#endif // SYCCQE3DFITS_SYSTEMATICS_H

