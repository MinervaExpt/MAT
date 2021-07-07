//File: Model.h
//Brief: A Model is an internally consistent combination of several Reweighters.
//       Model's entire user interface is to be fed to Universe::GetWeight().
//       Internally, Model reorganizes how Reweighters are used to optimize performance
//       where possible.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_MODEL_H
#define PLOTUTILS_MODEL_H

//reweighters includes
#include "PlotUtils/reweighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Model
  {
    using weighters_t = std::vector<std::unique_ptr<Reweighter<UNIVERSE, EVENT>>>;

    public:
      Model(weighters_t&& weighters)
      {
        //Check for incompatible Reweighters
        for(auto whichWeighter = weighters.begin(); whichWeighter != weighters.end(); ++whichWeighter)
        {
          for(auto otherWeighter = std::next(whichWeighter); otherWeighter != weighters.end(); ++otherWeighter)
          {
            if(!(*whichWeighter)->IsCompatible(**otherWeighter) || !(*otherWeighter)->IsCompatible(**whichWeighter)) throw std::runtime_error("Reweighters named " + (*whichWeighter)->GetName() + " and " + (*otherWeighter)->GetName() + " are not compatible!");
          }
        }

        //Sort Reweighters by whether they need to check IsVerticalOnly()
        for(auto& weighter: weighters)
        {
          if(weighter->DependsReco()) fDependsRecoWeighters.push_back(std::move(weighter));
          else fDependsTruthOnlyWeighters.push_back(std::move(weighter));
        }
      }

      void SetEntry(const UNIVERSE& univ, const EVENT& event)
      {
        fCVTruthOnlyWeight = 1;
        for(const auto& weighter: fDependsTruthOnlyWeighters) fCVTruthOnlyWeight *= weighter->GetWeight(univ, event);

        fCVDependsRecoWeight = 1;
        for(const auto& weighter: fDependsRecoWeighters) fCVDependsRecoWeight *= weighter->GetWeight(univ, event);
      }

      double GetWeight(const UNIVERSE& univ, const EVENT& event) const
      {
        double weight = fCVTruthOnlyWeight; //TODO: From event?

        if(!univ.IsVerticalOnly())
        {
          for(const auto& weighter: fDependsRecoWeighters) weight *= weighter->GetWeight(univ, event);
        }
        else weight *= fCVDependsRecoWeight;

        //If we ever have systematics that shift truth quantities,
        //do the same kind of IsVerticalOnly() check for truth here.

        return weight * univ.GetWeightRatioToCV();
      }

    private:
      weighters_t fDependsRecoWeighters;
      weighters_t fDependsTruthOnlyWeighters;

      //TODO: Put these in EVENT instead?  Otherwise, I need a separate function
      //      to set them once per event.
      double fCVTruthOnlyWeight;
      double fCVDependsRecoWeight;
  };
}

#endif //PLOTUTILS_MODEL_H
