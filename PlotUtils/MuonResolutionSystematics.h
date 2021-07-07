#ifndef MuonResolutionSYSTEMATICS_H
#define MuonResolutionSYSTEMATICS_H

#include "TreeWrapper.h"
#include "NSFDefaults.h"
#include "PlotUtilsPhysicalConstants.h"

/* instantiate this uncertainty with
 
 SystMap muonR_systematics = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain,NSFDefaults::muonResolution_Err);
  error_bands.insert(muonR_systematics.begin(), muonR_systematics.end());
 
 */

namespace PlotUtils{
  template<class T>
  class MuonResolutionUniverse: public T
  {
    public:
      MuonResolutionUniverse(typename T::config_t chw, double nsigma , double fracUncertainty = NSFDefaults::muonResolution_Err );

      double GetMuonResolutionMomentumShift() const;
      virtual double GetPmu() const /*override*/;

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_fracUncertainty;
  };

  template<class T>
  class MuonAngleXResolutionUniverse: public T
    {
    public:
      MuonAngleXResolutionUniverse(typename T::config_t chw, double nsigma , double fracUncertainty = NSFDefaults::muon_angle_frac_res);
      
      double GetTrueThetaXmu() const /*override*/;
      double GetThetaXmuShift() const /*override*/;
      virtual double GetThetaXmu() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_fracUncertainty;
    };

    template<class T>
    class MuonAngleYResolutionUniverse: public T
    {
    public:
      MuonAngleYResolutionUniverse(typename T::config_t chw, double nsigma , double fracUncertainty = NSFDefaults::muon_angle_frac_res);
      
      double GetTrueThetaYmu() const /*override*/;
      double GetThetaYmuShift() const /*override*/;
      virtual double GetThetaYmu() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_fracUncertainty;
    };
}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "MuonResolutionSystematics.cxx"

#endif // MuonResolutionSYSTEMATICS_H

