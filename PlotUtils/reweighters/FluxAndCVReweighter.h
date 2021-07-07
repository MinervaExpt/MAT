//File: FluxAndCVReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_FluxAndCVREWEIGHTER_H
#define PLOTUTILS_FluxAndCVREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/FluxReweighter.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class FluxAndCVReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      FluxAndCVReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~FluxAndCVReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        int nu_pdg = 14;
        if(univ.IsPlaylistME(univ.GetPlaylist())) nu_pdg = univ.GetInt("mc_incoming");

        return PlotUtils::flux_reweighter(univ.GetPlaylist(), nu_pdg, univ.UseNuEConstraint(),
                                    univ.GetNFluxUniverses()).GetFluxCVWeight(univ.GetDouble("mc_incomingE") * 1e-3, nu_pdg);
      }

      std::string GetName() const override { return "FluxAndCV"; }
      bool DependsReco() const override { return false; }
  };
}

#endif //PLOTUTILS_FluxAndCVREWEIGHTER_H
