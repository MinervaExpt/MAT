//File: GeantNeutronCVReweighter.h
//Brief: A Reweighter that changes the inelastic interaction cross section for FS neutrons.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_GeantNeutronCVREWEIGHTER_H
#define PLOTUTILS_GeantNeutronCVREWEIGHTER_H

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/GeantHadronSystematics.h"

//Reweighter includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class GeantNeutronCVReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      GeantNeutronCVReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~GeantNeutronCVReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        if(univ.IsTruth()) return 1; // No efficiency reweighting for truth events

        univ.SetupMHRWeighter(); //TODO: do this once in the constructor
        return PlotUtils::weight_hadron<PlotUtils::TreeWrapper*>().reweightNeutronCV(univ);
      }

      std::string GetName() const override { return "GeantNeutronCV"; }
      bool DependsReco() const override { return false; }
  };
}

#endif //PLOTUTILS_GeantNeutronCVREWEIGHTER_H
