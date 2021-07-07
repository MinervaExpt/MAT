//File: MKReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_MKREWEIGHTER_H
#define PLOTUTILS_MKREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/weightMK.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MKReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      MKReweighter(): Reweighter<UNIVERSE, EVENT>(), fCalculator(std::string(std::getenv("MPARAMFILESROOT")) + "/data/Reweight/output_ratio_genie_neut_for_MKmodel.root")
      {
      }

      virtual ~MKReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        return fCalculator.getWeight(univ.GetTree(), univ.GetEntry());
      }

      std::string GetName() const override { return "MK"; }
      bool DependsReco() const override { return false; }

      virtual bool IsCompatible(const Reweighter<UNIVERSE, EVENT>& other) const
      {
        //TODO: Is the MK reweight compatible with our RPA reweight?  Low Q2 pion suppression?  Coherent reweight?
        return other.GetName() != "LowRecoil2p2hTune";
      }

    private:
      mutable PlotUtils::weightMK fCalculator; //TODO: Why can't GetWeight() be const here?  I don't see any obvious reasons at first glance, but I'm far from an expert on this.
  };
}

#endif //PLOTUTILS_MKREWEIGHTER_H
