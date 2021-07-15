//File: Cutter.cxx
//Brief: A Cutter keeps track of the counters needed to report efficiency and purity
//       after an event selection.  It keeps track of 4 sets of Cuts: reco pre-cuts, reco sideband cuts, truth signal
//       definition, and truth phase space.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//N.B.: Including the header in the .cxx AND the .cxx in the header
//      can be solved more efficiently than with these include guards.
//      It should be solved in the build system instead.
//      However, cmt is not smart enough to deal with that.
#ifndef PLOTUTILS_CUTTER_CXX
#define PLOTUTILS_CUTTER_CXX

//PlotUtils includes
#include "Cutter.h"
#include "Table.h"

//c++ includes
#include <algorithm>
#include <cassert>

namespace MAT
{
  template <class UNIVERSE, class EVENT>
  Cutter<UNIVERSE, EVENT>::Cutter(reco_t&& recoPre, reco_t&& recoSideband,
                 truth_t&& truthSignal, truth_t&& truthPhaseSpace): fSumAllAnaToolWeights(0), fSumSignalTruthWeights(0),
                                                                    fSumSignalAnaToolWeights(0), fSumAllTruthWeights(0), fNRecoEntries(0), fNTruthEntries(0), 
                                                                    fRecoPreCuts(std::move(recoPre)), fRecoSidebandCuts(std::move(recoSideband)),
                                                                    fTruthSignalDef(std::move(truthSignal)), fTruthPhaseSpace(std::move(truthPhaseSpace)),
                                                                    fLastEntryCVChecked(-1)
  {
  }

  template <class UNIVERSE, class EVENT>
  std::bitset<64> Cutter<UNIVERSE, EVENT>::isMCSelected(const UNIVERSE& univ, EVENT& event, double weight)
  {
    bool isSignalForCuts = false;
    getMCStats(univ, weight, isSignalForCuts, isCV(univ));
    if(univ.IsVerticalOnly())
    {
      if((univ.GetEntry() != fLastEntryCVChecked) || isCV(univ)) //I must always call checkSelection() on the CV to get a useful cut table.
                                                                 //You will call checkSelection() on isVertical() universes at most twice
                                                                 //with this optimization.
      {
        fCVPassedCuts = checkSelection(univ, event, weight, isSignalForCuts);
        fLastEntryCVChecked = univ.GetEntry();
      }
      return fCVPassedCuts;
    }

    return checkSelection(univ, event, weight, isSignalForCuts);
  }

  template <class UNIVERSE, class EVENT>
  std::bitset<64> Cutter<UNIVERSE, EVENT>::isDataSelected(const UNIVERSE& univ, EVENT& event)
  {
    if(isCV(univ)) ++fNRecoEntries; //TODO: This is not true if I do any cuts before calling isSelected().
                                    //      Think of Heidi's column-wise cuts.
    return checkSelection(univ, event);
  }

  template <class UNIVERSE, class EVENT>
  std::bitset<64> Cutter<UNIVERSE, EVENT>::checkSelection(const UNIVERSE& univ, EVENT& event, double weight, const bool isSignalForCuts)
  {
    if(!std::all_of(fRecoPreCuts.begin(), fRecoPreCuts.end(), [&univ, &event, weight, isSignalForCuts](std::unique_ptr<Cut<UNIVERSE, EVENT> >& cut)
                                                              { return cut->passesCut(univ, event, weight, isSignalForCuts); }))
    {
      return 0; //All sideband cuts set to false
    }

    constexpr int nCutsMax = 64;
    assert(1 + fRecoSidebandCuts.size() < nCutsMax && "Cutter needs to be upgraded to support more than 63 sideband cuts.");
    std::bitset<nCutsMax> result;
    result.set(); //Set all Cuts to passed by default so user can just check isSelected().all()

    int whichCut = 0;
    for(auto& cut: fRecoSidebandCuts)
    {
      if(!cut->passesCut(univ, event, weight, isSignalForCuts))
      {
        result.set(whichCut, false);
        weight = 0; //Keep checking for sideband, but this event has not been selected.
      }
      ++whichCut;
    }

    return result;
  }

  template <class UNIVERSE, class EVENT>
  void Cutter<UNIVERSE, EVENT>::getMCStats(const UNIVERSE& univ, double& weight, bool& isSignalForCuts, const bool isThisTheCV)
  {
    //If univ is the CV, keep statistics for summarize()
    if(isThisTheCV)
    {
      ++fNRecoEntries; //TODO: This is not true if I do any cuts before calling isSelected().
                       //      Think of Heidi's column-wise cuts.
                                                                                              
      fSumAllAnaToolWeights += weight;
                                                                                              
      isSignalForCuts = isSignal(univ) && isPhaseSpace(univ);
      if(isSignalForCuts) fSumSignalAnaToolWeights += weight;
    }
    else
    {
      weight = 0;
    }
  }

  template <class UNIVERSE, class EVENT>
  void Cutter<UNIVERSE, EVENT>::getTruthStats(const UNIVERSE& univ, double& weight, bool& isSigPS, const bool isThisTheCV)
  {
    //If univ is the CV, keep statistics for summarize()
    if(isThisTheCV)
    {
      ++fNTruthEntries; //TODO: This is not true if I do any cuts before calling isSelected().
                       //      Think of Heidi's column-wise cuts.
      
      fSumAllTruthWeights += weight;
      //Pass the weight to isSignal and isPhaseSpace to do the counting
      isSigPS = isSignal(univ,weight) && isPhaseSpace(univ,weight);
      if(isSigPS) fSumSignalTruthWeights += weight;
    }
    else
    {
      weight = 0;
    }
  }

  template <class UNIVERSE, class EVENT>
  bool Cutter<UNIVERSE, EVENT>::isSignal(const UNIVERSE& univ,const double weight)
  {
    return std::all_of(fTruthSignalDef.begin(), fTruthSignalDef.end(), [&univ,weight](const std::unique_ptr<SignalConstraint<UNIVERSE> >& def) { return def->passes(univ,weight); });
  }

  template <class UNIVERSE, class EVENT>
  bool Cutter<UNIVERSE, EVENT>::isPhaseSpace(const UNIVERSE& univ, const double weight)
  {
    return std::all_of(fTruthPhaseSpace.begin(), fTruthPhaseSpace.end(), [&univ,weight](const std::unique_ptr<SignalConstraint<UNIVERSE> >& def) { return def->passes(univ,weight); });
  }

  template <class UNIVERSE, class EVENT>
  bool Cutter<UNIVERSE, EVENT>::isEfficiencyDenom(const UNIVERSE& univ, double weight)
  {
    bool isSigPS = false;
    getTruthStats(univ, weight, isSigPS, isCV(univ));
    return isSignal(univ) && isPhaseSpace(univ);
  }

  template <class UNIVERSE, class EVENT>
  std::ostream& Cutter<UNIVERSE, EVENT>::summarizeReco(std::ostream& printTo) const
  {
    assert(fNRecoEntries > 0 && "No entries got into the cut table!  Perhaps you don't have a universe with ShortName() == \"cv\"?");

    //2 kinds of summaries are available: with Truth tree and without.
    //If a loop over the Truth tree has been run, then I can report efficiency
    //and purity.
    //Otherwise, I report just fraction of all events passed.  The latter is
    //more useful than it seems: it can be compared to data.
    if(fSumSignalTruthWeights > 0)
    {
      //Summary with Truth tree
      util::Table<6> truthSummary({"Cut Name",
                                   "Events",
                                   "\% Eff",
                                   "\% Purity",
                                   "Relative \% Eff",
                                   "Relative \% All"});
      truthSummary.appendRow("AnaTool",
                             fSumAllAnaToolWeights,
                             fSumSignalAnaToolWeights / fSumSignalTruthWeights * 100.,
                             fSumSignalAnaToolWeights / fSumAllAnaToolWeights * 100.,
                             fSumSignalAnaToolWeights / fSumSignalTruthWeights * 100.,
                             fSumAllAnaToolWeights / fSumAllTruthWeights * 100.);

      double prevSignal = fSumSignalAnaToolWeights;
      double prevTotal = fSumAllAnaToolWeights;

      auto summarizeCut = [&truthSummary, &prevSignal, &prevTotal, this](const std::unique_ptr<Cut<UNIVERSE, EVENT> >& cut)
                          {
                            truthSummary.appendRow(cut->getName(),
                                                   cut->getTotalPassed(),
                                                   cut->getSignalPassed() / this->fSumSignalTruthWeights * 100.,
                                                   cut->getSignalPassed() / cut->getTotalPassed() * 100.,
                                                   cut->getSignalPassed() / prevSignal * 100.,
                                                   (double)cut->getTotalPassed() / prevTotal * 100.);
                            prevSignal = cut->getSignalPassed();
                            prevTotal = cut->getTotalPassed();
                          };

      std::for_each(fRecoPreCuts.begin(), fRecoPreCuts.end(), summarizeCut);
      std::for_each(fRecoSidebandCuts.begin(), fRecoSidebandCuts.end(), summarizeCut);

      return truthSummary.print(printTo);
    }

    //else implied by return in above if()
    //Summary with only reco tree
    util::Table<3> recoSummary({"Cut Name",
                                "Total Entries",
                                "Relative \% Sample Left"});

    recoSummary.appendRow("AnaTool",
                          fNRecoEntries,
                          100.);

    int prevSampleLeft = fNRecoEntries;
    auto summarizeCut = [&recoSummary, &prevSampleLeft](const std::unique_ptr<Cut<UNIVERSE, EVENT> >& cut)
                        {
                          recoSummary.appendRow(cut->getName(),
                                                cut->getTotalPassed(),
                                                (double)cut->getTotalPassed() / prevSampleLeft * 100.);
                          prevSampleLeft = cut->getTotalPassed();
                        };

    std::for_each(fRecoPreCuts.begin(), fRecoPreCuts.end(), summarizeCut);
    std::for_each(fRecoSidebandCuts.begin(), fRecoSidebandCuts.end(), summarizeCut);

    return recoSummary.print(printTo);
  }

  template <class UNIVERSE, class EVENT>
  std::ostream& Cutter<UNIVERSE, EVENT>::summarizeTruth(std::ostream& printTo) const
  {
    util::Table<1> signalDef({"Signal Constraint"});
    for(const auto& sig: fTruthSignalDef) signalDef.appendRow(sig->getName());
    signalDef.print(printTo) << "\n\n";

    util::Table<1> phaseSpace({"Phase Space Constraint"});
    for(const auto& phase: fTruthPhaseSpace) phaseSpace.appendRow(phase->getName());
    return phaseSpace.print(printTo);
  }

  template <class UNIVERSE, class EVENT>
  std::ostream& Cutter<UNIVERSE, EVENT>::summarizeTruthWithStats(std::ostream& printTo) const
  {
    util::Table<3> signalDef({"   Signal Constraint   ",
                              "Total Entries",
                              "Relative \% Sample Left"});
    signalDef.appendRow("No Constraint",
                          fNTruthEntries,
                          100.);
    double prevSampleLeftSig = fNTruthEntries;
    auto summarizeCutSig = [&signalDef, &prevSampleLeftSig](const std::unique_ptr<SignalConstraint<UNIVERSE> >& constraints)
                        {  
                          signalDef.appendRow(constraints->getName(),
                                                constraints->getTotalPassed(),
                                                (double)constraints->getTotalPassed() / prevSampleLeftSig * 100.);
                          prevSampleLeftSig = constraints->getTotalPassed();
                        };

    std::for_each(fTruthSignalDef.begin(), fTruthSignalDef.end(), summarizeCutSig);

    signalDef.print(printTo) << "\n\n";


    util::Table<3> phaseSpace({"Phase Space Constraint",
                              "Total Entries",
                              "Relative \% Sample Left"});
    phaseSpace.appendRow("After Signal Constraint",
                          prevSampleLeftSig,
                          100.);
    double prevSampleLeftPS = prevSampleLeftSig;
    auto summarizeCutPS = [&phaseSpace, &prevSampleLeftPS](const std::unique_ptr<SignalConstraint<UNIVERSE> >& constraints)
                        {  
                          phaseSpace.appendRow(constraints->getName(),
                                                constraints->getTotalPassed(),
                                                (double)constraints->getTotalPassed() / prevSampleLeftPS * 100.);
                          prevSampleLeftPS = constraints->getTotalPassed();
                        };
    
    std::for_each(fTruthPhaseSpace.begin(), fTruthPhaseSpace.end(), summarizeCutPS);

    return phaseSpace.print(printTo);
  }

  template <class UNIVERSE, class EVENT>
  double Cutter<UNIVERSE, EVENT>::totalWeightPassed() const
  {
    if(!fRecoSidebandCuts.empty()) return fRecoSidebandCuts.back()->getTotalPassed();
    else if(!fRecoPreCuts.empty()) return fRecoPreCuts.back()->getTotalPassed();
    else return fSumAllAnaToolWeights;
  }

  template <class UNIVERSE, class EVENT>
  std::ostream& operator <<(std::ostream& printTo, const Cutter<UNIVERSE, EVENT>& printMe)
  {
    printMe.summarizeReco(printTo) << "\n\n";
    return printMe.summarizeTruth(printTo) << "\n\n";
  }

  template <class UNIVERSE, class EVENT>
  std::ostream& operator <(std::ostream& printTo, const Cutter<UNIVERSE, EVENT>& printMe)
  { 
    return printMe.summarizeTruthWithStats(printTo) << "\n\n";
  }

  template <class UNIVERSE, class EVENT>
  void Cutter<UNIVERSE, EVENT>::resetStats()
  {
    for(auto& cut: fRecoPreCuts) cut->resetStats();
    for(auto& cut: fRecoSidebandCuts) cut->resetStats();

    fNRecoEntries = 0;
    fNTruthEntries = 0;
    fSumAllAnaToolWeights = 0;
    fSumSignalAnaToolWeights = 0;
    fSumAllTruthWeights = 0;
    fSumSignalTruthWeights = 0;
  }

  template <class UNIVERSE, class EVENT>
  bool Cutter<UNIVERSE, EVENT>::isCV(const UNIVERSE& univ)
  {
    //TODO: Don't do a string comparison because this is usually called in a tight loop.
    //      I could imagine either taking a pointer to the CV or making this a DCVU function.
    //      For the DCVU function, I'd write a "systematic universe" that overrides nothing
    //      except IsCV().
    return !strcmp(univ.ShortName().c_str(), "cv");
  }


  //Expert interface
  //If you know whether a particular universe is the CV or have a way to find out
  //that's better than isCV(), call this function.
  template <class UNIVERSE, class EVENT>
  std::bitset<64> Cutter<UNIVERSE, EVENT>::isMCSelectedCV(const UNIVERSE& univ, EVENT& event, double weight)
  {
    assert(isCV(univ) && "This is an expert interface.  You must always pass the CV to it.");
                                                                                                                    
    bool isSignalForCuts = false;
    getMCStats(univ, weight, isSignalForCuts, true);
                                                                                                                    
    fCVPassedCuts = checkSelection(univ, event, weight, isSignalForCuts);
    return fCVPassedCuts;
  }

  //Expert interface                                                                                                
  //Don't keep any stats because I know I already got the CV.
  //Also, check Cuts for a group of UNIVERSES that all pass the same Cuts.
  template <class UNIVERSE, class EVENT>
  std::bitset<64> Cutter<UNIVERSE, EVENT>::isSelectedWithNoStats(const std::vector<UNIVERSE*>& univs, EVENT& event)
  {
    return checkSelection(*univs.front(), event, 0, false);
  }
}

#endif //PLOTUTILS_CUTTER_CXX
