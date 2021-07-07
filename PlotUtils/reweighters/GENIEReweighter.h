//File: GENIEReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_GENIEREWEIGHTER_H
#define PLOTUTILS_GENIEREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/GenieSystematics.cxx" //IsNonResPi()
#include "PlotUtils/MnvTuneSystematics.cxx" //IsCCRes()

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class GENIEReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      GENIEReweighter(bool useNonResPi, bool useDeuteriumPionTune): Reweighter<UNIVERSE, EVENT>(), fUseNonResPiReweight(useNonResPi), fUseDeuteriumPionTune(useDeuteriumPionTune)
      {
        //TODO: Deprecate and remove these from PlotUtils because GetRequiredUniverses() does the same job but with a lot less of a mess.
        UNIVERSE::SetNonResPiReweight(useNonResPi);
        UNIVERSE::SetDeuteriumGeniePiTune(useDeuteriumPionTune);
      }

      virtual ~GENIEReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        double nonResPiWgt = fUseNonResPiReweight && PlotUtils::IsNonResPi(univ)
                           ? PlotUtils::kNonResPiWeight : 1.;
        double deutWgt = fUseDeuteriumPionTune && PlotUtils::IsCCRes(univ) ?
                         ( PlotUtils::GetGenieParReweight(univ,"truth_genie_wgt_MaRES",
                                                        NSFDefaults::DEUTERIUM_MaRES,
                                                        NSFDefaults::GENIE_MaRES,
                                                        NSFDefaults::GENIE_MaRES_1Sig ) * 
                          NSFDefaults::DEUTERIUM_RES_NORM ) : 1.;
        return nonResPiWgt * deutWgt;
      }

      std::string GetName() const override { return "GENIE"; }
      bool DependsReco() const override { return false; }

    private:
      bool fUseNonResPiReweight;
      bool fUseDeuteriumPionTune;
  };
}

#endif //PLOTUTILS_GENIEREWEIGHTER_H
