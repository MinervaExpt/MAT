//File: LowQ2PiReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_LowQ2PiREWEIGHTER_H
#define PLOTUTILS_LowQ2PiREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/MnvTuneSystematics.cxx" //IsCCRes()

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class LowQ2PiReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      LowQ2PiReweighter(const std::string& channel): Reweighter<UNIVERSE, EVENT>(), fChannel(channel)
      {
      }

      virtual ~LowQ2PiReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        //variation 0 is the CV
        if(!PlotUtils::IsCCRes(univ)) return 1;
        return PlotUtils::weight_lowq2pi().getWeight(univ.GetQ2True() * 1e-6 /*GeV^2*/,
                                                     fChannel, 0);
      }

      std::string GetName() const override { return "LowQ2Pi_" + fChannel; }
      bool DependsReco() const override { return false; }

      std::vector<UNIVERSE*> GetRequiredUniverses() const override
      {
        return std::vector<UNIVERSE*>{}; //TODO: Return the low Q2 pion universes here with the right channel
      }

    private:
      std::string fChannel; //Which low Q2 pion reweight are we using?
  };
}

#endif //PLOTUTILS_LowQ2PiREWEIGHTER_H
