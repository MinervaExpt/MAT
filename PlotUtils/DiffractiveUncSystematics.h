#ifndef DIFFRACTIVEUNCSYSTEMATICS_H
#define DIFFRACTIVEUNCSYSTEMATICS_H

#include "TreeWrapper.h"
#include "NSFDefaults.h"

namespace PlotUtils{
  //Right now, just hardcoding this.  
  //To estimate uncertainty on diffractive pion events, putting an uncertainty on coherent events
  //  to simulate our lack of knowledge about diffractive events 
  //  50% for tracker, 20% for nuclear targets
  
  static double fracDiffractiveTrackerUnc = 0.5;
  static double fracDiffractiveTargetUnc  = 0.2;

  template<class T>
  class DiffractiveUncUniverse: public T
  {
    public:
      DiffractiveUncUniverse(typename T::config_t chw, double nsigma, double fracTrkUnc = fracDiffractiveTrackerUnc , double fracTarUnc = fracDiffractiveTargetUnc );
      
      virtual double GetDiffractiveUncWeight() const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return true; }/*override*/
    private:
      double m_fracTrkUnc;
      double m_fracTarUnc;
  };
  
}

// Explicit Instantiation, needed for loading into python
//! Make sure you put this into Ana/PlotUtils/dict/PlotUtilsDict.h
//template class

// For template classes, the header needs to know about the definitions
#include "DiffractiveUncSystematics.cxx"

#endif // DIFFRACTIVEUNCSYSTEMATICS_H

