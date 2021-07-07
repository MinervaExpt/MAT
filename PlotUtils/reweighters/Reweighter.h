//File: Reweighter.h
//Brief: A Reweighter changes the CV model into a different model using just a multiplicative
//       constant.  All vertical systematics are implemented by taking ratios to such weights.
//       Some Reweighters are mutually exclusive, and others are only needed for specific systematics.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_REWEIGHTER_H
#define PLOTUTILS_REWEIGHTER_H

#if __cplusplus < 201103L
  #define override
#endif

//c++ includes
#include <string>
#include <vector>

namespace PlotUtils
{
  namespace detail
  {
    struct empty;
  }

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Reweighter
  {
    public:
      Reweighter() = default;
      virtual ~Reweighter() = default;

      virtual double GetWeight(const UNIVERSE& univ, const EVENT& event) const = 0;
      virtual std::string GetName() const = 0;

      virtual bool DependsReco() const = 0;
      /*virtual bool DependsTruth() const;*/ //Not needed as of time of writing.

      virtual bool IsCompatible(const Reweighter& /*other*/) const { return true; }
      virtual std::vector<UNIVERSE*> GetRequiredUniverses() const { return std::vector<UNIVERSE*>{}; }
  };
}

#endif //PLOTUTILS_REWEIGHTER_H
