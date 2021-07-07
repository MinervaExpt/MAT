//File: LowRecoil2p2hReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_LowRecoil2p2hREWEIGHTER_H
#define PLOTUTILS_LowRecoil2p2hREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/MnvTuneSystematics.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class LowRecoil2p2hReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      LowRecoil2p2hReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~LowRecoil2p2hReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        //variation 0 is the CV
        return PlotUtils::GetLowRecoil2p2hWeight(univ, univ.Getq0True() / 1000 /* GeV */,
                                           univ.Getq3True() / 1000 /* GeV */,
                                           0);
      }

      std::string GetName() const override { return "LowRecoil2p2hTune"; }
      bool DependsReco() const override { return false; }
  };
}

#endif //PLOTUTILS_LowRecoil2p2hREWEIGHTER_H
