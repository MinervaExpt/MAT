#ifndef BASECVUNIVERSE_H
#define BASECVUNIVERSE_H

#include <iostream>
#include "ChainWrapper.h"
#include "FluxReweighter.h"
#include "Math/Vector4D.h"
#include "NSFDefaults.h"
#include "PlotUtilsPhysicalConstants.h"


namespace PlotUtils {
  class BaseUniverse {
    protected:
      // MOVE THESE TO A PHYSICSCALCULATORS.H HEADER
      //! helper functions for calculating physics angles from measured angles
      double phi3D(double thetaX, double thetaY) const;
      double theta3D(double thetaX, double thetaY) const;
      double calcq0(const double Enu, const double Elep) const;
      double calcq3(const double Q2, const double Enu, const double Elep) const;

      // member variables
      PlotUtils::TreeWrapper* m_chw; // chain of events
      double m_nsigma;               // n sigma that universe event shifts
      Long64_t m_entry;              // current chain entry
      static bool m_is_truth;
      
    public:
      BaseUniverse() {}
  
      typedef PlotUtils::ChainWrapper* config_t;
  
      //! Constructor
      BaseUniverse(PlotUtils::TreeWrapper* config, double nsigma = 0);
  
      //! Destructor
      virtual ~BaseUniverse() {}

      //===============================================
      // Qualities of an abstract Universe
      //===============================================
      virtual bool IsVerticalOnly() const { return false; }

      //Hook for systematic universes that explicitly modify weights to override
      virtual double GetWeightRatioToCV() const { return 1; }

    protected:
      // Set per-entry member variables here.  This function
      // is called in SetEntry(), so you are guaranteed not
      // to be stuck with "stale" event information.
      virtual void OnNewEntry() {}

    public:
      // Some systematics (e.g. muon efficiency) should only be applied to reco,
      // hence the need for m_is_truth
      static void SetTruth(bool is_truth) { m_is_truth = is_truth; }
      bool IsTruth() const { return m_is_truth; }
  
      void SetEntry(Long64_t entry) {
        m_entry = entry;
        OnNewEntry();
      }
      int GetEntry() const { return m_entry; }
      double GetSigma() const { return m_nsigma; }
      TreeWrapper* GetTree() const { return m_chw; } //TODO: I don't understand how this can be const.  It's needed for code to use this with the MK reweighter.
  
      virtual std::string ShortName() const { return "cv"; }
      virtual std::string LatexName() const { return "Central value"; }
      virtual std::string GetAnaToolName() const {
        return m_chw->GetTree()->GetName();
      }

      //===============================================
      // Generic Branch Getters
      //===============================================
      //! Get Dummy Values
      double GetDummyVar() const { return -999.; }
      double GetDummyVar1Arg(const int idx1) const {
        (void)idx1;
        return -999.;
      }
      double GetDummyVar2Arg(const int idx1, const int idx2) const {
        (void)idx1;
        (void)idx2;
        return -999.;
      }
      int GetInt(const char* name) const;
      bool GetBool(const char* name) const;
      double GetDouble(const char* name) const;
      int GetVecElemInt(const char* name, const int i) const;
      double GetVecElem(const char* name, const int i) const;
      double GetVecElem(const char* name, const int i, const int j) const;
      std::vector<int> GetVecInt(const char* name) const;
      std::vector<double> GetVecDouble(const char* name) const;
      std::vector<std::vector<double> > GetVecOfVecDouble(const char* name) const;
      template <typename T>
      std::vector<T> GetVec(const char* name) const {
        return m_chw->GetValueVector<T>(name, m_entry);
      }
  };

}

#endif  // BASECVUNIVERSE
