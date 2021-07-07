//==============================================================================
//In this file several inclusive cuts are defined.
//
//Each cut is a class that inherits from PlotUtils::Cut.
//
//At minimum, cuts must have a name and they must override the checkCut
//function.
//
//You cut also takes two template parameters:
//* UNIVERSE template parameter: Cuts are built on Universe objects. Plug
//in your CVUniverse as a template parameter so the cut can access your
//branches and do so correctly within a systematic universe.
//
// * EVENT template parameter: sometimes cuts need to do more than just return
//a bool when you call the checkCut function. The event object is a
//user-specifiable object that can hold onto information that is learned within
//checkCut. E.g. A HasMichel cut determines which tracks are potential pion
//tracks. See CCPionCuts.h for a fleshed out example making use of EVENT.
//==============================================================================

//Original Author: Andrew Olivier aolivier@ur.rochester.edu

#include "PlotUtils/Cut.h"
#include <sstream>
#ifndef BEN_CCINCLUSIVECUTS_H
#define BEN_CCINCLUSIVECUTS_H

#ifdef __GCCXML__
#define override
namespace std {
  std::string to_string(double x) {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }
}
#endif

namespace reco
{
  //============================================================================
  //Example 1: The simplest cut example. Just derive from Cut<> base class.
  //============================================================================
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasMINOSMatch: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      // Constructor
      HasMINOSMatch(): PlotUtils::Cut<UNIVERSE, EVENT>("Has MINOS Match") {}

    private:
      // THE cut function
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        // Call a CVUniverse member function to make the cut
        return univ.IsMinosMatchMuon();
      }
  };

  //============================================================================
  //Example 2: Many cuts require variables to be above, below, or at some value.
  //To simply help avoid typos and reduce a little typing, use Minimum,
  //Maximum, and IsSame helper templates.
  //In this case the muon energy is required to be at least X, at minimum.
  //The specific value of X gets set when you instantiate this cut (in
  //getCCInclusiveCuts).
  //============================================================================
#ifndef __GCCXML__
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  using MuonEnergyMin = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetEmu, EVENT>;
#endif
  //============================================================================
  //Example 3: The first, commented-out dead time cut inherits directly from Ben's
  //CVUniverse class. (In form, it resembles case 1, by the way.) While this is
  //valid, it means that no one else can use this cut (because they don't also
  //use my CVUniverse).
  //============================================================================
  /*class NoDeadtime: public Cut<CVUniverse, detail::empty>
  {
    public:
      NoDeadtime(const int nDead = 1): Cut("Dead Discriminators"), m_nDeadAllowed(nDead)
      {
      }

    private:
      const int m_nDeadAllowed; //Number of dead discriminators allowed

      bool checkCut(const CVUniverse& univ, const detail::empty&) const override
      {
        return univ.GetTDead() < m_nDeadAllowed;
      }
  };*/

  //============================================================================
  //Example 3 contd: It's better to make the cut a template so that any analyzer's
  //CVUniverse will work with it.
  //Furthermore, look how much space we save by using the Maximum helper
  //template.
  //============================================================================
#ifndef __GCCXML__
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  using NoDeadtime = PlotUtils::Maximum<UNIVERSE, int, &UNIVERSE::GetTDead, EVENT>;
#endif
  //============================================================================
  //Example 4: Another case like Example 1.
  //Physics-wise, we're making sure that the event is from a neutrino as
  //opposed to an antineutrino.  Negative muon charge means neutrino
  //interaction.
  //============================================================================
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsNeutrino: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      IsNeutrino(): PlotUtils::Cut<UNIVERSE, EVENT>("Muon Charge Sign") {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetMuonQP() < 0; 
      }
  };

  //Other Cuts I need to reproduce Dan's ME 2D inclusive analysis
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MaxMuonAngle: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      MaxMuonAngle(const double max): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Muon Theta < ") + std::to_string(max)), fMax(max*M_PI/180.) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetThetamu() < fMax;
      }

      const double fMax; //radians
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ZRange: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetVertexZ() >= fMin && univ.GetVertexZ() <= fMax;
      }

      const double fMin;
      const double fMax;
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Apothem: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    Apothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem), fSlope(-1./sqrt(3.))//A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        const ROOT::Math::XYZTVector vertex = univ.GetVertex();
        return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
               && (fabs(vertex.x()) < fApothem);
      }

      const double fApothem;
      const double fSlope; 
  };
#ifndef __GCCXML__
  //============================================================================
  // This function instantiates each of the above cuts and adds them to a
  // container, over which we'll loop during our event selection to apply the
  // cuts.
  // 
  // The return type for this function is a `cuts_t<UNIVERSE, EVENT>`, which is
  // shorthand for std::vector<std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT>>>;
  //============================================================================
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  PlotUtils::cuts_t<UNIVERSE, EVENT> GetCCInclusiveCuts()
  {
    PlotUtils::cuts_t<UNIVERSE, EVENT> inclusive_cuts;
    inclusive_cuts.emplace_back(new HasMINOSMatch<UNIVERSE, EVENT>());
    inclusive_cuts.emplace_back(new MuonEnergyMin<UNIVERSE, EVENT>(2e3, "Emu"));
    inclusive_cuts.emplace_back(new NoDeadtime<UNIVERSE, EVENT>(1, "Deadtime"));
    inclusive_cuts.emplace_back(new IsNeutrino<UNIVERSE, EVENT>());

    return inclusive_cuts;
  }

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  PlotUtils::cuts_t<UNIVERSE, EVENT> GetCCInclusive2DCuts()
  { 
    PlotUtils::cuts_t<UNIVERSE, EVENT> inclusive_cuts;
    inclusive_cuts.emplace_back(new ZRange<UNIVERSE, EVENT>("Tracker", 5980, 8422));
    inclusive_cuts.emplace_back(new Apothem<UNIVERSE, EVENT>(850.));
    inclusive_cuts.emplace_back(new MaxMuonAngle<UNIVERSE, EVENT>(20.));
    inclusive_cuts.emplace_back(new HasMINOSMatch<UNIVERSE, EVENT>());
    inclusive_cuts.emplace_back(new NoDeadtime<UNIVERSE, EVENT>(1, "Deadtime"));
    inclusive_cuts.emplace_back(new IsNeutrino<UNIVERSE, EVENT>());

    return inclusive_cuts;
  }
#endif
  //TODO: MnvGENIEv1, nu-e constraint, 50 flux universes
  //TODO: Binning: [ 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60]
  //               [ 0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5]

}
#ifdef __GCCXML__
#undef override
#endif

#endif //BEN_CCINCLUSIVECUTS_H
