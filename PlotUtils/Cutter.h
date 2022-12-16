//File: Cutter.h
//Brief: A Cutter keeps track of the counters needed to report efficiency and purity
//       after an event selection.  It keeps track of 4 sets of Cuts: reco pre-cuts, reco sideband cuts, truth signal
//       definition, and truth phase space.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_CUTTER_H
#define PLOTUTILS_CUTTER_H

//PlotUtils includes
#include "PlotUtils/Cut.h"

//c++ includes
#include <iostream>
#include <memory>
#include <bitset>

namespace evt
{
  class CVUniverse;
}

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Cutter
  {
    public:
      using reco_t = std::vector<std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT> > >;
      using truth_t = std::vector<std::unique_ptr<PlotUtils::SignalConstraint<UNIVERSE> > >;

      Cutter(reco_t&& recoPre, reco_t&& recoSideband, truth_t&& truthSignal, truth_t&& truthPhaseSpace);

      //Hooks to keep statistics within your event loop.

      //Look for isSelected().all() for the main analysis.
      //isSelected()[n] is the status of the nth sideband cut.
      //isSelected().none() means the pre-cuts failed and no sidebands were checked to save computation time.
      std::bitset<64> isMCSelected(const UNIVERSE& univ, EVENT& event, double weight); //Extra statistics for MC
      std::bitset<64> isDataSelected(const UNIVERSE& univ, EVENT& event); //For data where I can't keep statistics

      //isSignal() && isPhaseSpace() with some extra bookkeeping to report overall efficiency and purity.
      bool isEfficiencyDenom(const UNIVERSE& univ, double weight);

      //No statistics kept for these
      bool isSignal(const UNIVERSE& univ, const double weight = 0);
      bool isPhaseSpace(const UNIVERSE& univ, const double weight = 0);

      std::ostream& summarizeReco(std::ostream& printTo) const;
      std::ostream& summarizeTruth(std::ostream& printTo) const;
      std::ostream& summarizeTruthWithStats(std::ostream& printTo) const;
      double totalWeightPassed() const; //Get total weight that passed all cuts

      void resetStats();

      //Extra expert interface.  Assumes that all univs pass the same Cuts.
      //std::bitset<64> isMCSelected(const std::vector<UNIVERSE*>& univs, EVENT& event, double weight);
      std::bitset<64> isMCSelectedCV(const UNIVERSE& univ, EVENT& event, double weight);
                                                                                                                  
      //Extra expert interface.  Check Cuts with no statistics kept when I know that I've already checked the CV.
      std::bitset<64> isSelectedWithNoStats(const std::vector<UNIVERSE*>& univs, EVENT& event);

    private:
      double fSumAllAnaToolWeights;
      double fSumSignalTruthWeights;
      double fSumSignalAnaToolWeights;
      double fSumAllTruthWeights;
      long int fNRecoEntries;
      long int fNTruthEntries;

      reco_t fRecoPreCuts;
      reco_t fRecoSidebandCuts;
      truth_t fTruthSignalDef;
      truth_t fTruthPhaseSpace;

      //Have I checked whether the CV passes all cuts for this event?  If so, I can assume the same
      //for all IsVertical() universes.
      std::bitset<64> fCVPassedCuts;
      EVENT fCVEvent;
      int fLastEntryCVChecked;
      std::bitset<64> checkSelection(const UNIVERSE& univ, EVENT& event, double weight = 1., const bool isSignalForCuts = true);
      void getMCStats(const UNIVERSE& univ, double& weight, bool& isSignalForCuts, const bool isThisTheCV);
      void getTruthStats(const UNIVERSE& univ, double& weight, bool& isSig, const bool isThisTheCV);
      bool isCV(const UNIVERSE& univ); //Hook for a more efficient CV check later
  };

  template <class UNIVERSE, class EVENT>
  std::ostream& operator <<(std::ostream& printTo, const Cutter<UNIVERSE, EVENT>& printMe);
  template <class UNIVERSE, class EVENT>
  std::ostream& operator <(std::ostream& printTo, const Cutter<UNIVERSE, EVENT>& printMe);
}

#include "Cutter.cxx"

#endif //PLOTUTILS_CUTTER_H
