//File: FSIReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_FSIREWEIGHTER_H
#define PLOTUTILS_FSIREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/weight_fsi.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

//c++ includes
#include <cassert>

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class FSIReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      FSIReweighter(bool useElastic, bool useAbsorption): Reweighter<UNIVERSE, EVENT>(), fUseElastic(useElastic), fUseAbsorption(useAbsorption)
      {
        assert((useElastic || useAbsorption) && "Attempted to use the GENIE FSI bug fix reweight with neither elastic nor absorption reweighting!");

        fCalculator.UseTrackingThreshold();
      }

      virtual ~FSIReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        fCalculator.calcWeights(univ.GetInt("mc_incoming"), univ.GetInt("mc_primaryLepton"), univ.GetInt("mc_charm"), univ.GetInt("mc_intType"), univ.GetInt("mc_targetA"), univ.GetInt("mc_targetZ"), univ.GetInt("mc_resID"), univ.GetInt("mc_er_nPart"), univ.GetVecInt("mc_er_ID"), univ.GetVecInt("mc_er_status"), univ.GetVecInt("mc_er_FD"), univ.GetVecInt("mc_er_LD"), univ.GetVecInt("mc_er_mother"), univ.GetVecDouble("mc_er_Px"), univ.GetVecDouble("mc_er_Py"), univ.GetVecDouble("mc_er_Pz"), univ.GetVecDouble("mc_er_E"), univ.GetEntry());

        double weight = 1.;
        if(fUseElastic) weight *= fCalculator.GetElasticWeight(1);
        if(fUseAbsorption) weight *= fCalculator.GetAbsorptionWeight();
        return weight;
      }

      std::string GetName() const override { return "GENIE_FSI_Bug_Fix"; }
      bool DependsReco() const override { return false; }

    private:
      bool fUseElastic;
      bool fUseAbsorption;

      mutable PlotUtils::weight_fsi fCalculator; //TODO: I'll bet GetElasticWeight() and GetAbsorptionWeight() could both be const
  };
}

#endif //PLOTUTILS_FSIREWEIGHTER_H
