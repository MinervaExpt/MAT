//File: RPAReweighter.h
//Brief: A Reweighter that reweights events based on FS kinematics to simulate the random phase approximation for nucleon initial state kinematics.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_RPAREWEIGHTER_H
#define PLOTUTILS_RPAREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/MnvTuneSystematics.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class RPAReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      RPAReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~RPAReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        //variation 0 is the CV
        return PlotUtils::GetRPAWeight(univ, univ.Getq0True() / 1000 /* GeV */,
                                 univ.Getq3True() / 1000 /* GeV */, 0,
                                 univ.IsProcessingNX());
      }

      std::string GetName() const override { return "RPA"; }
      bool DependsReco() const override { return false; }
  };
}

#endif //PLOTUTILS_RPAREWEIGHTER_H
