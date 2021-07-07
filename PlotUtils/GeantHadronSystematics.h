#ifndef GeanHadronSystematics_h
#define GeanHadronSystematics_h

#include "MnvHadronReweight.h"
#include "NSFDefaults.h"
#include <cstdlib>

namespace PlotUtils {
/*! @brief Implements Geant uncertainty on hadron (pi, p, n) (re)interaction
           probabilities. MnvHadronReweight class does all physics hard work.

           MnvHadronReweight needs to renormalize in order to fix the event
           rate. There are a set of standard renormalization weights in 
           $PLOTUTILSROOT/data/mhrwKineRenorm, but the user may want to include
           different hadrons or different analysis volumes, which will require
           new renormalizations.   
*/
  template <class T>
  class GeantHadronUniverse : public T {
   public:
    GeantHadronUniverse(typename T::config_t chw, double nsigma, int pdg );

    virtual double GetGeantHadronWeight() const; /* override */
    double GetWeightRatioToCV() const;

    int m_pdg;
    std::string m_part_name;
    virtual bool IsVerticalOnly() const { return true; } /* override */
    virtual std::string ShortName() const {
      return "GEANT_" + m_part_name;
    } /* override */
    virtual std::string LatexName() const {
      return "GEANT " + m_part_name;
    } /* override */ 
  };
}  // namespace PlotUtils

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "GeantHadronSystematics.cxx"

#endif  // GeanHadronSystematics_h
