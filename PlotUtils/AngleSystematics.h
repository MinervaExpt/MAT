#ifndef AngleSYSTEMATICS_H
#define AngleSYSTEMATICS_H

#include "TreeWrapper.h"
#include "NSFDefaults.h"

namespace PlotUtils{
  ///////////////////////////////////////////////////////////////
  // Some notes: this is the result of a low nu study
  // which looked at the difference between the residual of MC and Data
  // which found a correction and error to the muon angle. So this is the error 
  // on the beam and and the error on the muon angle resolution

  template<class T>
  class BeamAngleXUniverse: public T
    {
    public:
      BeamAngleXUniverse(typename T::config_t chw, double nsigma , double Uncertainty = NSFDefaults::beamThetaX_Err);
      
      
      virtual double GetBeamAngleOffsetX() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_Uncertainty;
    };

    template<class T>
    class BeamAngleYUniverse: public T
    {
    public:
      BeamAngleYUniverse(typename T::config_t chw, double nsigma , double Uncertainty = NSFDefaults::beamThetaY_Err);
      
      
      virtual double GetBeamAngleOffsetY() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/
    private:
      double m_Uncertainty;
    };

  
}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "AngleSystematics.cxx"

#endif // AngleSYSTEMATICS_H

