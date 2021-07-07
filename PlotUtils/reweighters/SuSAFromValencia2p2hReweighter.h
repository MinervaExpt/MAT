//File: SuSAFromValencia2p2hReweighter.h
//Brief: A Reweighter that converts the Valencia 2p2h
//       model in our default Monte Carlo simulation
//       to the SuSA group's 2p2h prediction.  Note
//       that there is a separate weight (not yet
//       implemented) for the QE model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_SuSAFromValencia2p2hREWEIGHTER_H
#define PLOTUTILS_SuSAFromValencia2p2hREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/weightSusaValenciaClass.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class SuSAFromValencia2p2hReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      SuSAFromValencia2p2hReweighter(): Reweighter<UNIVERSE, EVENT>(),
                                        fWeighter(std::string(std::getenv("MPARAMFILESROOT"))+"/data/Reweight/outRatioSusaValencia2p2h-nu12C-9GeV-oP0-20191217.root")
      {
      }

      virtual ~SuSAFromValencia2p2hReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        if(univ.GetInt("mc_intType") == 8) return fWeighter.getWeight(univ.Getq0True() * 1.e-3, univ.Getq3True() * 1.e-3, 0);
        return 1;
      }

      std::string GetName() const override { return "SuSA2p2h"; }
      bool DependsReco() const override { return false; }

      //You must turn off the low recoil 2p2h tune (Phil tune) to reweight Valencia 2p2h to SuSA 2p2h.
      virtual bool IsCompatible(const Reweighter<UNIVERSE, EVENT>& other) const { return other.GetName() != "LowRecoil2p2hTune"; }

    private:
      weightSusaValenciaClass fWeighter;
  };
}

#endif //PLOTUTILS_SuSAFromValencia2p2hREWEIGHTER_H
