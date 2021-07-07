//File: AMUDISReweighter.h
//Brief: A Reweighter that modifies the DIS event rate using the AMU weighting.
//Comment: This could be adapted to allow for user input of which DIS weight
//         to use (see PlotUtils/weightDIS.h for the options).
//Author: David Last dlast@sas.upenn.edu

#ifndef PLOTUTILS_AMUDISREWEIGHTER_H
#define PLOTUTILS_AMUDISREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/weightDIS.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class AMUDISReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
    
    PlotUtils::weightDIS* m_weighter;

    AMUDISReweighter(): Reweighter<UNIVERSE, EVENT>()
    {
      m_weighter = new PlotUtils::weightDIS(PlotUtils::weightDIS::kAMU, true);
    }
    
    virtual ~AMUDISReweighter() = default;

    double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
    {
      double ret=1.0;
      const double Enu=univ.GetEnuTrue()/1000.0;//This is TruthFunction defined in MeV, but only needed in GeV.
	
      //Defined to return MeV^2 from the analysis CVUniverse true quantities of neutrino energy, muon energy, and muon angle
      //These may be slightly different than Dan's use in CCQENuInclusive
      const double Q2=univ.GetTrueExperimentersQ2();
      //Defined to return MeV^2 from the analysis CVUniverse true quantities of neutrino energy, muon energy, muon angle, and the above Q2
      //These may be slightly different than Dan's use in CCQENuInclusive
      const double W=univ.GetTrueExperimentersW();	
      
      //GetInteractionType() currently defined in our own CVUniverse, I imagine it'd be fine to move into SystCalcs/TruthFunctions.h at some point
      bool trueDIS=( (W/1000.0 >= 2.0) && (Q2/1.0e6 >= 1.0) && (univ.GetInteractionType() == 3) );
      if(!trueDIS) return 1.0; //Only reweight true DIS.
      
      //Bjorken True variables are also defined in our own CVUniverse, and would move easily to TruthFunctions.
      const double x=univ.GetBjorkenXTrue();
      const double y=univ.GetBjorkenYTrue();
      const int targetZ = univ.GetTargetZTrue(); //TruthFunctions defined      
      
      //PlotUtils::weightDIS* m_weighter = new PlotUtils::weightDIS(PlotUtils::weightDIS::kAMU);
      ret=m_weighter->getWeight(x,y,Enu,targetZ);
      //delete m_weighter;

      return ret;
    }

    std::string GetName() const override { return "AMUDIS"; }
    bool DependsReco() const override { return false; }
    /*
  private:
  PlotUtils::weightDIS* fWeighter;*/
      //Should keep an eye out on whether this would be needed here. 
      //I presume not given the exclusive nature of the intended low 
      //Q2 Pi use case

      /*
      std::vector<UNIVERSE*> GetRequiredUniverses() const override
      {
        return std::vector<UNIVERSE*>{}; //TODO: Return the low Q2 pion universes here with the right channel
      }
    */
  };
}

#endif //PLOTUTILS_AMUDISREWEIGHTER_H
